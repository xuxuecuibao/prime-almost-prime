"""Replay the three conditional proofs and verify overlap of the ranges."""
from decimal import Decimal
from hashlib import sha256
from pathlib import Path
import subprocess
import sys

if not __debug__:
    raise SystemExit("Enable Python assertions.")

root = Path(__file__).resolve().parent

def replay(folder, program, args, expected_sha256):
    result = subprocess.run(
        [sys.executable, "-B", program, *args], cwd=root / folder,
        capture_output=True, text=True, check=True
    )
    lines = result.stdout.splitlines()
    # Match the reference line sequence independently of platform line endings.
    digest = sha256("\n".join(lines).encode("utf-8")).hexdigest()
    if digest != expected_sha256:
        raise ArithmeticError(f"Replay differs from the reference for {folder}/{program}")
    return lines

def finite(lines, indexes, expected_count):
    assert len(lines) == expected_count
    start = previous = None
    worst = Decimal(0)
    smallest = Decimal(1)
    for line in lines:
        cells = line.split()
        left, right = map(Decimal, cells[:2])
        bad = Decimal(cells[indexes[0]])
        mass = Decimal(cells[indexes[1]])
        assert right > left and (previous is None or left == previous)
        assert 0 <= bad < 1 and mass > 0
        if start is None:
            start = left
        previous = right
        worst = max(worst, bad)
        smallest = min(smallest, mass)
    print(start, previous, "rows", len(lines), "maximum displayed excluded ratio", worst,
          "minimum displayed candidate coefficient", smallest, flush=True)
    return start, previous

low = finite(replay("low", "low.py", ["167"],
    "27d810f49ba71830e2ebff8aee8c9c49096b51ccb4f1fdcb1f34b0b7b98c821b"), (4,5), 125)
middle = finite(replay("middle", "verify.py", ["165","900","5"],
    "56f3d33d8ac651231f2765468148d3d649212a0db3c4f820942312523ab2c4b1"), (2,3), 147)
assert low[0] == Decimal("42.75")
assert low[0].exp() < 4*10**18 and middle[0] <= low[1] and middle[1] == 900
print(*replay("tail", "dimension.py", [],
    "42dbc6d45f1cfe8d32aa8d6371c7491db45e5f355e97cac0420633765b2162b2"),sep="\n")
print(*replay("tail", "check.py", [],
    "298c228a028f0cf1a84a97f62ebc186145a1220c6c379d545769975c1a5b17d9"),sep="\n")
print("All conditional bounds and interval joins passed.")
