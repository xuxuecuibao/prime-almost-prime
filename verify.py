"""Reproduce all four parts of the article (standard Python 3.11+, no options)."""
import subprocess
import sys
from pathlib import Path

root = Path(__file__).resolve().parent
if not __debug__:
    raise SystemExit("Run with assertions enabled, without -O.")
for folder, entry in [("p4cert", "verify.py"), ("p3ucert", "weighted.py"),
                      ("p3gcert", "verify_all.py"), ("p2cert", "verify.py")]:
    print(f"\n{folder}/{entry}", flush=True)
    subprocess.run([sys.executable, "-B", entry], cwd=root / folder, check=True)
