import os
from pathlib import Path
from dotenv import load_dotenv

load_dotenv(Path(__file__).parent / ".env")

api_key = os.environ.get("ALPHAGENOME_API_KEY", "")
if not api_key or api_key == "your_api_key_here":
    print("ERROR: ALPHAGENOME_API_KEY is not set. Edit .env and add your real key.")
    exit(1)

print(f"API key loaded: {api_key[:6]}{'*' * (len(api_key) - 6)}")

from alphagenome.models import dna_client
from alphagenome.data import genome
from alphagenome.models.dna_client import OutputType

print("Initializing AlphaGenome client...")
model = dna_client.create(api_key)
print(f"Client type: {type(model)}")

# Dummy variant: rs7903146 (TCF7L2, chr10:112998590, T>C) — well-known T2D locus
# Interval must be 1,048,576 bp (1 Mb) centered on the variant
CHROM = "chr10"
POS_1BASED = 112998590  # 1-based (VCF-style)
POS_0BASED = POS_1BASED - 1  # Interval uses 0-based
HALF = 524288
interval = genome.Interval(chromosome=CHROM, start=POS_0BASED - HALF, end=POS_0BASED + HALF)
variant = genome.Variant(chromosome=CHROM, position=POS_1BASED, reference_bases="T", alternate_bases="C")

print(f"\nInterval: {interval}")
print(f"Variant:  {variant}")
print("Submitting to AlphaGenome (this may take ~30s)...")

result = model.predict_variant(
    interval,
    variant,
    requested_outputs=[OutputType.RNA_SEQ],
    ontology_terms=None,
)

print("\n--- Response shape ---")
print(f"Result type:       {type(result)}")
print(f"Result attributes: {[a for a in dir(result) if not a.startswith('_')]}")
if hasattr(result, "tracks"):
    print(f"tracks type:       {type(result.tracks)}")
    print(f"tracks shape:      {getattr(result.tracks, 'shape', result.tracks)}")
