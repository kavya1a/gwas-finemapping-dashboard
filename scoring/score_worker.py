"""Subprocess-backed AlphaGenome scorer with real timeout cancellation.

The original ThreadPoolExecutor + future.result(timeout=...) + shutdown(wait=False)
pattern does not cancel a timed-out API call — shutdown(wait=False) only blocks
new submissions; the running gRPC call stays alive on the threadpool, holding
API quota and a TCP connection until the server eventually responds. Over a
long batch this leaks N concurrent in-flight calls and can trip rate limits.

Here each call runs in a subprocess. On timeout the parent terminate()s the
subprocess (killing the gRPC call) and respawns a fresh worker. The model is
built once per worker, not once per call, so serial throughput is comparable
to the original threadpool pattern.
"""

from __future__ import annotations

import multiprocessing
import queue as _queue
import threading


def _worker_main(api_key, req_q, res_q):
    try:
        from alphagenome.models import dna_client
        from scoring.composite import score_single_variant
        model = dna_client.create(api_key)
    except Exception as exc:
        try:
            res_q.put({"error": f"worker_init: {exc}"})
        except Exception:
            pass
        return

    while True:
        try:
            payload = req_q.get()
        except (EOFError, OSError):
            return
        if payload is None:
            return
        variant_input, tissue_profile = payload
        try:
            res_q.put(score_single_variant(model, variant_input, tissue_profile))
        except Exception as exc:
            res_q.put({"error": str(exc)})


class ScoreWorker:
    """Single subprocess worker with terminate-on-timeout semantics.

    score(...) returns the raw score_single_variant result on success, or
    {"error": "api_timeout"} / {"error": "<str>"} on failure. The worker is
    transparently re-spawned after a timeout.
    """

    def __init__(self, api_key: str) -> None:
        self._api_key = api_key
        self._ctx = multiprocessing.get_context("spawn")
        self._proc = None
        self._req_q = None
        self._res_q = None
        self._spawn()

    def _spawn(self) -> None:
        self._req_q = self._ctx.Queue()
        self._res_q = self._ctx.Queue()
        self._proc = self._ctx.Process(
            target=_worker_main,
            args=(self._api_key, self._req_q, self._res_q),
            daemon=True,
        )
        self._proc.start()

    def _kill(self) -> None:
        if self._proc is None:
            return
        if self._proc.is_alive():
            self._proc.terminate()
            self._proc.join(timeout=5)
            if self._proc.is_alive():
                self._proc.kill()
                self._proc.join()

    def score(self, variant_input, tissue_profile, timeout: float) -> dict:
        try:
            self._req_q.put((variant_input, tissue_profile))
        except Exception as exc:
            self._kill()
            self._spawn()
            return {"error": f"queue_put: {exc}"}
        try:
            return self._res_q.get(timeout=timeout)
        except _queue.Empty:
            self._kill()
            self._spawn()
            return {"error": "api_timeout"}

    def close(self) -> None:
        try:
            if self._req_q is not None:
                self._req_q.put(None)
        except Exception:
            pass
        self._kill()


class ScoreWorkerPool:
    """Pool of N ScoreWorkers, safe for threadpool-style parallel dispatch.

    score(...) is thread-safe; callers check out a worker for the duration of
    one call and return it afterwards. A timed-out worker is killed and
    re-spawned in place — the same pool slot is reused, so the pool size
    remains constant.
    """

    def __init__(self, n_workers: int, api_key: str) -> None:
        self._workers = [ScoreWorker(api_key) for _ in range(n_workers)]
        self._available: _queue.Queue = _queue.Queue()
        for w in self._workers:
            self._available.put(w)
        self._lock = threading.Lock()

    def score(self, variant_input, tissue_profile, timeout: float) -> dict:
        worker = self._available.get()
        try:
            return worker.score(variant_input, tissue_profile, timeout)
        finally:
            self._available.put(worker)

    def close(self) -> None:
        with self._lock:
            for w in self._workers:
                w.close()
            self._workers = []
