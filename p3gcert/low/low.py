"""GRH, uniform slabs: certified fixed-N prime plus a squarefree P3."""
from functools import lru_cache
from bisect import bisect_right
import math,sys
from directed import I,c,log,exp,mn,mx,primes
if not __debug__: raise RuntimeError('Run with assertions enabled')
def floor(x):
    assert int(x.lo)==int(x.hi)
    return int(x.lo)
P=primes(floor(exp(13)))
SCALE=10**28
prefix1=[0]; prefix2=[0]
Qs=[I(0)]; Ds=[I(1)]
for p in P:
    prefix1.append(prefix1[-1]+(SCALE+p-2)//(p-1))
    prefix2.append(prefix2[-1]+(SCALE+p*(p-1)-1)//(p*(p-1)))
    if len(Qs)<=150:
        Qs.append(Qs[-1]+log(p)); Ds.append(Ds[-1]*(1-1/I(p-1)))
@lru_cache(None)
def position(a):
    return 0 if a==0 else bisect_right(P,floor(exp(I(a))))
@lru_cache(None)
def RS(a,b):
    a,b=I(a),I(b)
    # sum_{p>exp(a)}1/(p(p-1)) <= 1/floor(exp(a)).
    return log(b/a)+1/(2*a**2)+1/(2*b**2)+1/(exp(a)-1)
def checked(lo,hi,n,millih):
    l=I(lo);u=I(hi);h=I(millih)/1000
    det=4*h-1; mu=I('.006')/(l-3)
    assert (l-I('2.45')).lo > log(599).hi
    assert (2*(I('.04')*u**2+u)*exp(-l/2)/c).hi < mu.lo
    f=exp(-l/2);f2=exp(-l)
    err0=(2+(n+2)*u)*f2/c
    line=I('2.55')*f/c
    def error(v):
        return line*(v/2+2)+err0
    s=1-mu-2**n*error(Qs[n]); assert s.lo>0
    # The best branch is selected by ordinary floats; its validity and value
    # are independently checked by interval arithmetic before use.
    qs=[float(q.hi) for q in Qs[:n+1]]
    ratios=[d/Ds[n] for d in Ds[:n+1]]
    rs=[float(r.hi) for r in ratios]; cif=float(c.lo)
    linef=float(line.hi);errf=float(err0.hi);mu0=float(mu.hi)
    def cap(v):
        vv=float(v.hi); ev=math.exp(vv);best=(1e300,0,'grh')
        for i,(q,r) in enumerate(zip(qs,rs)):
            a=r*(1+mu0+2**i*ev*(linef*((q+vv)/2+2)+errf))
            if a<best[0]:best=(a,i,'grh')
            if q+vv<lo:
                b=2*r/(cif*(1-(q+vv)/lo))
                if b<best[0]:best=(b,i,'bt')
        _,i,mode=best
        if mode=='grh':return ratios[i]*(1+mu+2**i*exp(v)*error(Qs[i]+v))
        den=1-(Qs[i]+v)/l;assert den.lo>0
        return 2*ratios[i]/(c*den)
    marked=I(0);square=I(0)
    for j in range(13*8):
        a=j/8;b=(j+1)/8
        left=max(n,position(a));right=position(b)
        if right<=left:continue
        mass=I(prefix1[right]-prefix1[left])/SCALE
        mass2=I(prefix2[right]-prefix2[left])/SCALE
        if (h-I(a)/u).hi>0:
            marked+=mx(0,h-I(a)/u)*cap(I(b))*mass
        square+=cap(I(2*b))*mass2
    a=13
    while (h*u).hi>a:
        b=a+1
        marked+=mx(0,h-I(a)/u)*cap(I(b))*RS(a,b)
        a=b
    # Above exp(13) squares: integers m range from 1 to <N, so no +sqrt(N) term.
    square+=u/(c*Ds[n]*(exp(13)-1))
    other=(u*exp(-l/2)+2*u*exp(-l))/(c*Ds[n])
    bad=(marked/det+square+other)/s
    return bad.hi,s.lo
if __name__=='__main__':
    end=int(sys.argv[1]) if len(sys.argv)>1 else 154
    start=float(sys.argv[2]) if len(sys.argv)>2 else 42.75
    left=start
    while left<end:
        right=min(end,math.floor(left)+1)
        n=max(9,int(.49*left-10))
        h=400 if left<55 else (425 if left<89 else 450)
        B,S=checked(left,right,n,h)
        print(left,right,n,h,format(B,'.12f'),format(S,'.12f'),flush=True)
        assert B<1 and S>0
        left=right
