"""Tests for scoring/score_worker.py — no API calls, all via injected stubs.

The worker subprocess uses `spawn` context, so factories must be importable
module-level functions (pickle stores them by qualified name). Each factory
takes api_key and returns a callable (variant_input, tissue_profile) -> dict.
"""

from __future__ import annotations

import concurrent.futures
import time

import pytest

from scoring.score_worker import ScoreWorker, ScoreWorkerPool


# ---------------------------------------------------------------------------
# Stub factories (must be module-level for pickling under spawn)
# ---------------------------------------------------------------------------

def _echo_factory(api_key):
    def score(vi, profile):
        return {"rsid": vi, "composite_score": 0.5, "error": None}
    return score


def _slow_factory(api_key):
    def score(vi, profile):
        time.sleep(5)  # well past any test timeout
        return {"composite_score": 0.5}
    return score


def _crashing_factory(api_key):
    raise RuntimeError("simulated init failure")


def _raising_factory(api_key):
    def score(vi, profile):
        raise ValueError("simulated per-call failure")
    return score


# ---------------------------------------------------------------------------
# ScoreWorker
# ---------------------------------------------------------------------------

def test_success_returns_factory_payload():
    w = ScoreWorker(api_key="x", score_factory=_echo_factory)
    try:
        result = w.score("rs1", None, timeout=10.0)
        assert result == {"rsid": "rs1", "composite_score": 0.5, "error": None}
    finally:
        w.close()


def test_timeout_kills_subprocess_and_respawns():
    w = ScoreWorker(api_key="x", score_factory=_slow_factory)
    try:
        original_pid = w._proc.pid
        assert w._proc.is_alive()

        result = w.score("rs1", None, timeout=0.5)
        assert result == {"error": "api_timeout"}

        # Original subprocess should be reaped; a fresh one in its place.
        assert w._proc.pid != original_pid
        assert w._proc.is_alive()
    finally:
        w.close()


def test_close_reaps_subprocess():
    w = ScoreWorker(api_key="x", score_factory=_echo_factory)
    assert w._proc.is_alive()
    w.close()
    assert not w._proc.is_alive()


def test_worker_init_failure_surfaces_as_error():
    w = ScoreWorker(api_key="x", score_factory=_crashing_factory)
    try:
        result = w.score("rs1", None, timeout=3.0)
        assert "error" in result
        # Either the init-failure message reached the queue, or the subprocess
        # exited before posting and the parent saw a timeout — both are valid
        # surfacings of "worker did not produce a result".
        assert "worker_init" in result["error"] or result["error"] == "api_timeout"
    finally:
        w.close()


def test_per_call_exception_surfaces_as_error_dict():
    w = ScoreWorker(api_key="x", score_factory=_raising_factory)
    try:
        result = w.score("rs1", None, timeout=5.0)
        assert result == {"error": "simulated per-call failure"}
        # Worker stays alive after a per-call exception — only the call failed.
        assert w._proc.is_alive()
    finally:
        w.close()


def test_worker_serves_multiple_requests_on_one_subprocess():
    """Verify the model is built once and reused — same pid across calls."""
    w = ScoreWorker(api_key="x", score_factory=_echo_factory)
    try:
        pid_before = w._proc.pid
        for i in range(5):
            result = w.score(f"rs{i}", None, timeout=10.0)
            assert result["rsid"] == f"rs{i}"
        assert w._proc.pid == pid_before
    finally:
        w.close()


# ---------------------------------------------------------------------------
# ScoreWorkerPool
# ---------------------------------------------------------------------------

def test_pool_round_robin_under_concurrency():
    """N=3 workers handle 6 concurrent requests via a 3-thread dispatcher."""
    pool = ScoreWorkerPool(3, api_key="x", score_factory=_echo_factory)
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
            futures = [ex.submit(pool.score, f"rs{i}", None, 10.0)
                       for i in range(6)]
            results = [f.result() for f in futures]

        rsids = {r["rsid"] for r in results}
        assert rsids == {f"rs{i}" for i in range(6)}
    finally:
        pool.close()


def test_pool_timeout_on_one_worker_does_not_kill_others():
    """A timeout-induced kill+respawn affects only the offending worker;
    sibling workers in the pool continue serving requests normally."""
    pool = ScoreWorkerPool(2, api_key="x", score_factory=_echo_factory)
    try:
        # Capture all pids; one will change after we synthesize a timeout.
        original_pids = sorted(w._proc.pid for w in pool._workers)

        # Pull one worker directly and force a timeout on it.
        victim = pool._available.get()
        # Replace its factory by writing a slow result-handler is not possible
        # post-spawn — instead, simulate by sending a payload it can't answer
        # in time. Easiest: shrink the timeout to 0.01s so the echo doesn't
        # complete the round-trip in time.
        result = victim.score("rs-victim", None, timeout=0.001)
        assert result == {"error": "api_timeout"}
        pool._available.put(victim)

        # Remaining pool: 1 original-pid + 1 new (respawned) pid.
        new_pids = sorted(w._proc.pid for w in pool._workers)
        survivors = set(original_pids) & set(new_pids)
        assert len(survivors) == 1, (
            f"expected exactly 1 untouched worker, got pids "
            f"original={original_pids} now={new_pids}"
        )

        # Pool still services requests cleanly across both slots.
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
            futures = [ex.submit(pool.score, f"rs{i}", None, 10.0)
                       for i in range(4)]
            for f in futures:
                r = f.result()
                assert r.get("composite_score") == 0.5
    finally:
        pool.close()


def test_pool_close_reaps_all_workers():
    pool = ScoreWorkerPool(3, api_key="x", score_factory=_echo_factory)
    workers_snapshot = list(pool._workers)
    assert all(w._proc.is_alive() for w in workers_snapshot)
    pool.close()
    # Give terminate() a moment to settle on slower CI.
    time.sleep(0.2)
    assert all(not w._proc.is_alive() for w in workers_snapshot)
