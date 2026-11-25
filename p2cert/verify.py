"""Reproduce the GRH two-factor estimates in their order of use."""
from pathlib import Path
import subprocess
import sys

base = Path(__file__).resolve().parent
jobs = [
    ("triple", "verify.py", []),
    ("core", "verify.py", ["30000", ".485"]),
    ("core", "marked.py", []),
    ("shape", "certify.py", ["30000", ".10", ".485", "256", "96"]),
    ("core", "tail.py", []),
]
if not __debug__:
    raise SystemExit("Enable Python assertions.")
for folder, name, arguments in jobs:
    print(f"{folder}/{name}", flush=True)
    subprocess.run([sys.executable, "-B", name, *arguments],
                   cwd=base / folder, check=True)
