"""Verify the optional two large-log regimes in logarithmic arithmetic."""
from math import comb, factorial
from decimal import Decimal
from interval_math import I, log

if not __debug__:
    raise RuntimeError("Enable assertions")

theta = I(".485")
a = I(".1")
C = I(10**7)
lam = I("1.296")
for left, right, order in [(1300000, 10**8, 26), (10**8, None, None)]:
    L = I(left)
    if order is None:
        logM = 800000 * log(2)
        logq = 800000 * log(C)
        density = I(0)
        eta = I(0)
    else:
        logM = log(2**24 * sum(comb(800000,j) for j in range(order+1)))
        logq = 24 * log(97) + order * log(C)
        eta = lam**order / factorial(order)
        density = 50 * eta * I(right)**2
    assert (theta * L + logq).hi < L.lo
    assert density.hi < Decimal(".00001")
    # Work with logarithms of the positive exponential error terms.
    # Each decreases throughout the interval because (1/2-theta)*L>5.
    charerr = log(20) + logM + 29*log(10) + 5*log(L) - (I(".5")-theta)*L
    induced = log(800) + logM + 4*log(L) - a*L
    assert charerr.hi < -100 and induced.hi < -100
    # Relative mass errors for the base and switched sequence B_1:
    # (1.001 eta + M*20 L^2 exp(-aL) + M*10^29 L^4 exp(-L/2))/.03.
    anchor1 = logM + log(20) + 2*log(L) - a*L
    anchor2 = logM + 29*log(10) + 4*log(L) - L/2
    assert anchor1.hi < -100 and anchor2.hi < -100
    assert (I("1.001")*eta/I(".03")).hi < Decimal(".00001")
    # z^-1/(1-z^-1) <= 2 z^-1; the other three deletion terms are <= z^-1.
    deletions = log(100) + 2*log(L) - a*L
    assert deletions.hi < -100
    print("L >=",left,"upper",right,"J",order,"density at most",density.hi,
          "log character error <",charerr.hi,"log induced error <",induced.hi)
