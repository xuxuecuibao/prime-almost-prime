"""Check 18 adjoining four-factor intervals using the original rational sieve."""

from fractions import Fraction
from pathlib import Path
import sys

from mont_cert import cand, upper

rows = [
    (42, 45, 317, 48), (45, 49, 317, 49), (49, 53, 317, 49),
    (53, 57, 317, 50), (57, 61, 317, 51), (61, 65, 317, 52),
    (65, 68, 317, 53), (68, 70, 317, 54), (70, 71, 317, 54),
    (71, 72, 317, 54), (72, 73, 317, 55), (73, 78, 317, 56),
    (78, 83, 317, 57), (83, 87, 317, 58), (87, 89, 317, 58),
    (89, 90, 317, 59), (90, 91, 317, 59), (91, 92, 281, 59),
]

if not __debug__:
    raise RuntimeError("Run Python with assertions enabled.")

endpoint = 42
for left, right, prime_cutoff, hundred_h in rows:
    assert left == endpoint and left < right
    five_distinct, square = upper(left, right, prime_cutoff, hundred_h, 2)
    surplus = cand(left, prime_cutoff) - five_distinct - square
    assert surplus > Fraction(4811707, 10**9)
    print(left, right, (10**9 * surplus).__floor__())
    endpoint = right
assert endpoint == 92
