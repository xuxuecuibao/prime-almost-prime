"""Replay exactly the six families used for the unconditional four-factor result."""
import subprocess
import sys
from parameter_tables import SUPPORT_RECIPES, FINITE_DIFFERENCE_CHOICES, iter_parameter_rows
from pathlib import Path

if not __debug__:
    raise RuntimeError("Leave assertions enabled.")
root = Path(__file__).resolve().parent
from decimal import Decimal
from p4_fd_verify import I, run as finite_difference_run
from support_cert import check as support_run
from upper_iv import variable, bands, rows as band_rows
from four_variable import ROWS
from two_block import run as two_run, rows as two_rows

def contiguous_table(table, start, end, runner):
    last=start
    for row, choices in iter_parameter_rows(table):
        assert row[0]==last and row[1]>last
        last=row[1]
        runner(row, choices)
    assert last==end

subprocess.run([sys.executable,'-B',str(root/'check_four_low.py')],check=True)
contiguous_table(SUPPORT_RECIPES,91,115,support_run)
contiguous_table(FINITE_DIFFERENCE_CHOICES,115,212,finite_difference_run)
last=212
for row in ROWS:
    assert row[0]==last
    last=row[1]
    s,b=variable(row)
    print(row[0],row[1],'candidate >=',s,'bad <=',b,flush=True)
    assert s>0 and b<Decimal('.979725')
assert last==460
last=450
for left,right in band_rows():
    assert left==last
    last=right
    s,b=bands(left,right)
    print(left,right,'candidate >=',s,'bad <=',b,flush=True)
    assert s>0 and b<Decimal('.933814')
assert last==2800
last=2800
for row in two_rows():
    assert row[0]==last
    last=row[1]
    s,b=two_run(*row)
    print(row[0],last,'candidate >=',s,'bad <=',b,flush=True)
    assert s>0 and b<Decimal('.866405')
    if last>=23500: break
assert last>=23500
