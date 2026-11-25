"""Product bounds needed to apply Theorem 6 with an empty auxiliary set."""
from decimal import Decimal
from interval_math import I, exp, log, mx, primes

if not __debug__:
    raise RuntimeError("Run with assertions enabled.")

cutoff = 1500
eps = 1/I(300)
v = I(1)
for p in primes(cutoff):
    v *= 1 - 1/I(p-1)
cumulative = I(0)
best_start = None
excess = I(0)
for p in primes(7350):
    if p <= cutoff:
        continue
    loglog = log(log(p))
    start = loglog - cumulative
    best_start = start if best_start is None else mx(best_start, start)
    cumulative += log(I(p-1)/(p-2))
    excess = mx(excess, best_start + cumulative - loglog)
assert (exp(excess)-1).hi < Decimal('.002902')
assert exp(I('8.9')).hi < 7350
large = exp(I('.001614')+I('.000002964')/log(cutoff)+2/I(cutoff-2))-1
assert large.hi < eps.lo
assert (v*log(cutoff)/(1+eps)).lo > Decimal('.60')
print('finite relative excess < .002902; full excess < 1/300; density times log(z) > .60')
