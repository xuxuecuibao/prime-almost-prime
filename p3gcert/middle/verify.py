"""A GRH-only finite interval proof using elementary, short Bonferroni polynomials."""
import sys
from bisect import bisect_right
from decimal import ROUND_FLOOR,ROUND_CEILING,Decimal as D
from functools import lru_cache
from math import factorial,comb
from directed import I,c,exp,log,primes

if not __debug__: raise RuntimeError("Assertions must be enabled.")
# Every logarithmic mark boundary is an integer multiple of 1/MESH.
MESH=16
SCALE=10**28
M=10**6
PS=primes(M)
ARR=[0]
for p in PS: ARR.append(ARR[-1]+(SCALE+p-3)//(p-2))
MS=log(M)
@lru_cache(None)
def power_at(i): return exp(I(i)/MESH)
def floorhi(a):return int(a.hi.to_integral_value(rounding=ROUND_FLOOR))
def ceillo(a):return int(a.lo.to_integral_value(rounding=ROUND_CEILING))
def exact_floor(x):
    a=int(x.lo.to_integral_value(rounding=ROUND_FLOOR))
    assert a==floorhi(x)
    return a
@lru_cache(None)
def edge(i): return exact_floor(power_at(i))
def exact_mass(first,stop):
    return I(ARR[bisect_right(PS,stop)]-ARR[bisect_right(PS,first)])/SCALE
def analytic_mass(x,y):
    assert x.lo>2 and y.lo>=x.hi
    return (log(y/x)+(4*x+6)*exp(-x/2)/25)/(1-2*exp(-x))
@lru_cache(None)
def mass_between(a,b):
    # a,b encode log boundaries in mesh units, negative a means prime with index -a
    if a==b:return I(0)
    top=I(b)/MESH
    if a<0:
        first=PS[-a-1]; x=log(first)
    else:
        x=I(a)/MESH
        first=min(M,edge(a)) if x.lo<=MS.hi else M
    if top.hi<=MS.lo: return exact_mass(first,edge(b))
    if x.hi<MS.lo:return exact_mass(first,M)+analytic_mass(MS,top)
    return analytic_mass(x,top)
@lru_cache(None)
def term_count(k,j,b):
    # Rosser-Schoenfeld pi(x)<1.25506*x/log(x); pi(x)-1-k odd noncore primes.
    if (I(b)/MESH).hi < MS.lo:
        m=bisect_right(PS,edge(b))-k
    else: m=floorhi(I('1.26')*power_at(b)/(I(b)/MESH))-k
    assert m>=0
    return sum(comb(m,s) for s in range(min(m,j)+1)),m
@lru_cache(None)
def mark_mass(v):
    x=I(v)/MESH;y=I(v+1)/MESH
    a=power_at(v);b=power_at(v+1)
    if b.hi<=M:
        n=bisect_right(PS,edge(v+1))-bisect_right(PS,edge(v))
        return I(n)/(a-1)
    # Schoenfeld, Theorem 10: |theta(t)-t| <= sqrt(t)*log(t)^2/(8*pi), t>=599.
    assert a.lo>599
    return (b-a+(exp(y/2)*y**2+exp(x/2)*x**2)/25)/((a-1)*x)
def ceil_rat(j,k):return -(-j//k)
def previous_core_density(k):
    out=I(1)
    for p in PS[:k]:out*=1-1/I(p-1)
    return out
def check(l,u,k,bnum,hnum):
    assert 160<=l<u<=1000
    z=(u*bnum*MESH)//1000
    h=I(hnum)/1000
    assert z>0 and (I(z)/MESH).lo>log(PS[k-1]).hi and 4*h.lo>1
    pd=previous_core_density(k)
    lam=mass_between(-k,z); vi=exp(lam)
    inv_delta=vi/pd
    logW=k*log(PS[k-1])
    mu=I('.006')/(l-3)
    lhalf=exp(-I(l)/2);linv=exp(-l)
    assert (2*(I('.04')*u**2+u)*lhalf/c).hi<mu.lo
    def epsilon(q):
        return (I('2.55')*(q/2+2)*lhalf+2*linv+2*(q+1)*u*linv)/c
    def cutoff(j,v):
        mass=mass_between(-k,v)
        n,m=term_count(k,j,v)
        tr=I(0) if m<=j else exp(mass)*mass**(j+1)/factorial(j+1)
        return tr,n
    lows=[]
    for j in (1,3,5,7,9,11,13,15):
        tr,n=cutoff(j,z)
        if tr.hi>=1:continue
        lows.append( ((1-mu)*(1-tr)-2**k*vi*n*epsilon(logW+j*I(z)/MESH)).lo )
    s0=I(max(lows));assert s0.lo>0
    cand=[]
    # For marked q, upper = intercept + slope*exp(v). Bound all epsilon at q<=exp(h*u+1/MESH).
    last=ceil_rat(hnum*u*MESH,1000)
    for ki in sorted(set([0,3,6,k//2,k])):
        if ki>k:continue
        ratio=previous_core_density(ki)*inv_delta
        cand.append((((1+mu)*ratio/s0).hi,
                     (ratio*2**ki*epsilon(logW+I(last)/MESH)/s0).hi))
    start=ceillo(MESH*log(PS[k-1]))
    assert (I(start)/MESH).lo>log(PS[k-1]).hi
    for vp in set(start+(z-start)*i//12 for i in range(13)):
        ratio=exp(mass_between(vp,z))
        for j in (2,4,6,8,10,12):
            tr,n=cutoff(j,vp)
            cand.append(((ratio*(1+mu)*(1+tr)/s0).hi,
                         (vi*2**k*n*epsilon(logW+(j*vp+last)*I(1)/MESH)/s0).hi))
    # Pareto domination between conservative upper lines removes avoidable comparisons.
    live=[]
    for a,b in sorted(cand):
        if not live or b<live[-1][1]:live.append((a,b))
    total=I(0)
    for v in range(z,last):
        p=power_at(v+1).hi
        # Each endpoint operation in this shortcut is rounded up.
        cap=min((a+(b*p).next_plus()).next_plus() for a,b in live)
        weight=h-I(v)/(MESH*u)
        assert weight.lo>=0
        total+=mark_mass(v)*weight*I(cap)
    mark=total/(4*h-1)
    square=inv_delta/s0/c*(6/(power_at(z)-1)+u/(exp(I(l)/3)-1))
    deleted=inv_delta/s0/c*(u*lhalf+2*u*linv)
    return (mark+square+deleted).hi,s0.lo,mark.hi,square.hi

def params(l):
    if l<225:return 8,70,450
    if l<310:return 12,65,450
    if l<500:return 18,55,470
    if l<750:return 24,50,470
    return 24,45,470
if __name__=='__main__':
    l=int(sys.argv[1]) if len(sys.argv)>1 else 165
    end=int(sys.argv[2]) if len(sys.argv)>2 else 900
    stride=int(sys.argv[3]) if len(sys.argv)>3 else 5
    while l<end:
        u=min(end,l+stride)
        bad,s0,m,q=check(l,u,*params(l))
        unit=D('.000000001')
        print(l,u,bad.quantize(unit,rounding=ROUND_CEILING),s0.quantize(unit,rounding=ROUND_FLOOR),m.quantize(unit,rounding=ROUND_CEILING),q.quantize(unit,rounding=ROUND_CEILING),flush=True)
        assert bad<1
        l=u
