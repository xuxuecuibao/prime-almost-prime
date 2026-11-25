"""Uniform GRH estimates for the switched three-prime sequence for log N >= 10000."""
from directed import I, exp, log, sqrt
from decimal import Decimal

if not __debug__:
    raise RuntimeError("Enable assertions.")
L=I(10000)
delta=I(".001")
a=I(".09999"); b=I(".33332")
t=L/3-2*delta
assert ((L/10-2*delta)/(L+6*delta)).lo>a.hi
assert (t/(L+6*delta)).lo>b.hi
# All lower endpoints and maximum sum of the endpoint logs in a supported cell.
assert (L/10-2*delta).lo > (I(".099")*L).hi
assert (2*delta/(L/10-2*delta)).hi<Decimal(".0001")
# Relative error of the integral of one tent against theta, from Schoenfeld RH.
assert (2*(L+2*delta)**2*exp(2*delta)*exp(-L/20+delta)/(25*delta)).hi<Decimal(".0001")
# Maximal explicit-formula error for one normalized nonprincipal cell transform.
unit=(24/delta*(L/2+2)+2*delta+L+2*delta+L)/(L/10-2*delta)
assert unit.hi < 10**6
assert (I(1002)**3*10**18*exp(3*delta)).hi < Decimal('1e29')
assert (6*I('1.0001')**3*exp(8*delta)/log(2)).hi<20
# Decreasing, nonnegative integrand from a to 1-2b. Directed left rectangles.
end=1-2*b
n=10000
step=(end-a)/n
upper=I(0)
for j in range(n):
    x=a+j*step
    upper += step*2*log((1-x-b)/b)/(x*(1-x))
assert upper.hi < Decimal('.9825')
mass=I('1.0001')**6 * exp(6*delta)/(1-1/t) * upper
assert mass.hi<Decimal('.990')
print("Logarithmic integral upper bound:",upper.hi)
print("Weighted symmetric triple mass / (N/log N) less than:",mass.hi)
print("Error per reduced class at d <= N: 1e29 L^3 sqrt(N) + 20 L N^.9 / phi(d)")
