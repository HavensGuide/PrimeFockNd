"""Run smoke benchmarks from an editable checkout."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from causal_defect_engine.benchmarks import run_all

if __name__ == "__main__":
    results = run_all()
    for name, result in results.items():
        print(f"{name}: {result.seconds:.4f}s {result.metrics}")
