"""Outward checks for the marked presieved masses. Place beside core/verify.py."""
from decimal import Decimal
from math import comb, factorial
from interval_math import I, log, exp

for left, order in [(30000, 22), (1300000, 26), (10**8, None)]:
    ell = I(left)
    if order is None:
        log_m = 800000 * log(2)
        eta = I(0)
    else:
        log_m = log(2**24 * sum(comb(800000, i) for i in range(order + 1)))
        eta = I('1.296') ** order / factorial(order)
    marked = log_m + 29 * log(10) + 4 * log(ell) - ell/6
    principal = log_m + log(20) + 2 * log(ell) - I('.1') * ell
    assert marked.hi < -100 and principal.hi < -100
    bound = (I('1.001') * eta + 2 * exp(-100)) / I('.03')
    assert bound.hi < Decimal('.00001')
    print("L >=", left, "marked relative error <", bound.hi)
