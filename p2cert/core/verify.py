"""Sieve hypothesis, pre-sieve truncation and all error budgets in Chen's inequality (GRH)."""
import sys
from decimal import Decimal
from math import comb, factorial
from interval_math import I, log, exp, primes

if not __debug__:
    raise RuntimeError("Run without -O")
Lmin = I(sys.argv[1] if len(sys.argv)>1 else 30000)
theta = I(sys.argv[2] if len(sys.argv)>2 else '.485')
Lmax = I(1300000)
a = I('.1')
zlog = a*Lmin
C = I(10**7)
B = I(97)
J = 22
kc = len(primes(97))
assert kc == 24
assert (I('1.26')*C/log(C)).hi < 800000
M = 2**kc*sum(comb(800000,i) for i in range(J+1))
mass = log(log(C)/log(B)) + 1/(2*log(B)**2) + 1/(2*log(C)**2) + 1/(B-1)
eta = mass**J/factorial(J)
excess = exp(I('.4')/log(C)**3 + 2/(C-2))-1
assert excess.hi < Decimal('.0001')
assert (I('.5')/log(C)).lo > Decimal('.03')
assert (I(kc)*log(B)+J*log(C)).hi < ((1-theta)*Lmin).lo
assert zlog.lo>log(C).hi and (theta-I(1)/3).lo>a.hi
assert (.5-theta).lo>0 and ((.5-theta)*Lmin).lo>5
assert zlog.lo>4 and (I('.5')*Lmin).lo>5
# |class mean - nominal/phi| <= K for each reduced residue class; principal
# loss for the switched triple <= J0/phi; K0=K/sqrt(N), J00=J0/N.
K0 = I(10**29)*Lmin**3
j00 = 20*Lmin*exp(-a*Lmin)
# Anchoring the two pre-sieved mass estimates, in N/L units divided by V_C.
basis = (eta*I('1.001') + M*Lmin*(j00 + K0*exp(-Lmin/2)))/I('.03')
c = (1-exp(-20))*(1-exp(-I('.001')))/I('.001')
assert c.lo > Decimal('.9994') and c.hi<1
assert basis.hi < Decimal('.0001')
assert (c-basis).lo>Decimal('.998')
assert (c+basis).hi<Decimal('1.001')
assert (I('.989257')+basis).hi<Decimal('.990')
# All base/marked/switched sieve remainders divided by N V_N(z)/L.
resd = 50*eta*Lmax**2
resc = 20*M*I(10**29)*Lmin**5*exp(-(.5-theta)*Lmin)
resp = 800*M*Lmin**4*exp(-a*Lmin)
# m=1, a|N, prime powers; powers of marked primes; switched small primes and diagonal.
exception = 20*Lmin**2*(exp(-Lmin/2)+exp(-Lmin/3)+exp(-(1-theta)*Lmin)
                         +1/(exp(zlog)-1))
remainder = resd + resc + resp + exception
assert remainder.hi<Decimal('.001')
print('Lmin', Lmin.lo, 'Lmax', Lmax.lo, 'theta',theta.lo)
print('lambda, eta, logM upper', mass.hi, eta.hi, log(M).hi)
print('core exponential excess upper',excess.hi)
print('peak anchor loss (N V_C/L units) <',basis.hi)
print('density, character, induced terms (N V_N(z)/L units) <',resd.hi,resc.hi,resp.hi)
print('total error including deletions <',remainder.hi)
