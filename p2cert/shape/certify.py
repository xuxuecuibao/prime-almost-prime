"""Sieve main terms and the unrestricted sorted three-prime mass, N>=exp(L0)."""
from interval_math import I,log,exp,primes
from decimal import Decimal as D,ROUND_FLOOR,ROUND_CEILING
from functools import lru_cache
import sys
if not __debug__:raise RuntimeError("Run with assertions enabled")
LEFT=I(int(sys.argv[1]) if len(sys.argv)>1 else 30000)
a=I(sys.argv[2] if len(sys.argv)>2 else '.1')
T=I(sys.argv[3] if len(sys.argv)>3 else '.485')
b=1/I(3); eps=I('.0001'); C=I(106); C_lower=I(108)
G=2*exp(I('.5772156','.5772157'))
M=2048
Jlo=[I(0)];Jhi=[I(0)]
for k in range(2*M):
    for v,col in ((I(k)/M,Jlo),(I(k+1)/M,Jhi)):
        col.append(col[-1]+log(1+v)/(2+v)/M)
def ftu(s):
    assert 1<=s.lo<=s.hi<=5
    p=max(0,int(((s.hi-3)*M).next_plus().to_integral_value(rounding=ROUND_CEILING)))
    return G/s*(1+Jhi[p])
def ftlo(s):
    assert 4<=s.lo<=s.hi<=6
    # Integral I(u)/(3+u), 0<=u<=s-4, paid on complete grid cells.
    n=max(0,int(((s.lo-4)*M).next_minus().to_integral_value(rounding=ROUND_FLOOR)))
    extra=sum((Jlo[k]/(3+I(k+1)/M)/M for k in range(n)),I(0))
    return G/s*(log(s-1)+extra)
def h(x):
    if x.hi<=2:return exp(-2)
    if x.lo>=3:return 3*exp(-x)/x
    if x.lo>=2:return exp(-x)
    return exp(-2)
def err(x,c=C):return c*eps*exp(2)*h(x)
def har(x,y):
    # Restoring endpoints (closed intervals), pg−(p-1)^−1 adjustment, etc.
    return log(y/x)+(1/x**2+1/y**2)/(2*LEFT**2)+3/(exp(x*LEFT)-2)
def delratio(x,y):
    # Full-prime upper; excluding each p|N loses at most 1/(z−2) in log, at most 1/a such p.
    return (x/y)*exp((1/x**2+1/y**2)/(2*LEFT**2)+(3+1/a)/(exp(x*LEFT)-2))
assert (a*LEFT).lo>log(10**12).hi
K=int(sys.argv[4]) if len(sys.argv)>4 else 256
first=ftlo(T/a)-err(T/a,C_lower)
marked=I(0); jmass=I(0); sq=I(0)
for i in range(K):
    u=a+(b-a)*i/K; v=a+(b-a)*(i+1)/K
    ss=(T-v)/a
    # The F or h at each t>=ss are no larger throughout this branch. By positivity f<=1<=F, F'<=0.
    marked+=har(u,v)*(ftu(ss)+err(ss))/2
    top=(1-u)/2
    if top.hi<b.lo:continue
    # Sum in q2: log q2/L <=(1-u)/2. Account for both endpoints. q3<=x=N/(q1q2).
    # Use Rosser-Schoenfeld Thm 1 to upper bound pi(x); every actually counted box has x>=w.
    grid=int(sys.argv[5]) if len(sys.argv)>5 else 96
    inn=I(0)
    for j in range(grid):
        col=b+(top-b)*j/grid; end=b+(top-b)*(j+1)/grid
        inn+=har(col,end)/(1-v-end)
    jmass+=har(u,v)*inn
jmass*=1+I('1.5')/(b*LEFT)
# Unweighted n/log n for each sorted triple, then switched global upper sieve normalized by V_N(z).
sw=jmass*delratio(a,T)*(G+err(I(1)))/2
# S_A <=1.001 N/L for marked terms and >=.998 N/L initial.
sur=I('.998')*first-I('1.001')*marked-sw
for n,t in [('first',first),('marked',marked),('mass_sorted',jmass),('switched',sw),('surplus',sur)]:print(n,t.lo,t.hi)
if a.lo==D('.1') and T.lo==D('.485') and LEFT.lo>=30000 and K==256 and grid==96:
    assert (2*jmass).hi < D('.990')
    assert delratio(a,T).hi < (I('1.001')*a/T).lo
    independent_switch=I('.990')/4*I('1.001')*a/T*(G+err(I(1)))
    independent_sur=I('.998')*first-I('1.001')*marked-independent_switch
    print('independent_switch',independent_switch.lo,independent_switch.hi)
    print('independent_surplus',independent_sur.lo,independent_sur.hi)
    assert independent_sur.lo > D('.0087')
