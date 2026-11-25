"""Finite four-factor inequalities with Montgomery large sieve. Pure Python 3.11+."""
from fractions import Fraction as F
from math import prod,isqrt
from bisect import bisect_left
from functools import lru_cache
S=10**10
def ceil(a,b):return (a+b-1)//b
def up(x):return ceil(x.numerator,x.denominator)
def roundup(x): return F(up(S*x),S)
def logseries(z):
    # Upper bound for 2*atanh(z); 30 terms and a geometric tail.
    return (2*sum((z**(2*i+1)/(2*i+1) for i in range(30)),F(0))
            +2*z**61/(61*(1-z*z)))
LOG2=logseries(F(1,3))
@lru_cache(None)
def logupper(x):
    m=x.bit_length()-1; r=F(x,2**m)
    return roundup(m*LOG2+logseries((r-1)/(r+1)))
@lru_cache(None)
def expbounds(x):
    # Lower bound by truncation. A geometric bound for the tail is used
    # only for x <= 12 (the boundaries of the log grid).
    t=v=F(1)
    for j in range(1,121):
        t=t*x/j; v+=t
    upper=v+t*x/121/(1-x/122) if x<=12 else None
    return v,upper
def expfloor(x):
    v,w=expbounds(x)
    assert int(v)==int(w)
    return int(v)
n=200000
ispr=bytearray([1])*(n+1);ispr[:2]=bytes(2)
for p in range(2,isqrt(n)+1):
 if ispr[p]:ispr[p*p::p]=bytes((n-p*p)//p+1)
primes=[p for p in range(3,n+1,2) if ispr[p]]
vgrid=[expfloor(F(i,10)) for i in range(121)]
assert vgrid[-1]==162754
@lru_cache(None)
def const(B):
 P=[p for p in primes if p<=B]
 delta=prod(F(p-2,p-1) for p in P)
 H=prod(F(p+1,p) for p in P)
 U=roundup(sum(logupper(p)/(p-1) for p in P))
 wh=[0]*120
 for p in primes:
  if B<p<=vgrid[-1]:wh[bisect_left(vgrid,p)-1]+=ceil(S,p-1)
 return P,delta,H,U,wh
def cand(l,B):
 P,V,*_=const(B);vp=[F(1)]
 for p in P:vp.append(vp[-1]*F(p-2,p-1))
 def loss(d,ph,t,lower):
  err=10 if d>400000 else F(1,10000*ph)+F(1650*d,ph)/expbounds(F(l,2))[0]
  stop=vp[t]/ph if lower else (1-vp[t])/ph+err
  if err>=stop or d>400000:return stop
  total=err
  for j in range(t):
   total+=loss(d*P[j],ph*(P[j]-1),j,not lower)
   if total>=stop:return stop
  return total
 assert 2*(B+2)*l/(V*expbounds(F(l))[0]) < F(1,100000)
 return 1-loss(1,1,len(P),True)/V-F(1,100000)
def rational_logupper(x):
 return logupper(x.numerator)+logupper(1) - loglower(x.denominator)
def loglower(x):
 m=x.bit_length()-1;r=F(x,2**m);a=(r-1)/(r+1)
 series=lambda a:2*sum((a**(2*i+1)/(2*i+1) for i in range(30)),F(0))
 return F((S*(m*series(F(1,3))+series(a))).__floor__(),S)
@lru_cache(None)
def upper(l,u,B,h100,j,clo=F(67,100)):
 """Returns total bad in units c delta N, square and >=5 distinct contributions. 
 The nuisance numerator +1/(cYmin) is bounded by 1e-8, support entirely a>R."""
 h=F(h100,100);P,V,H,U,weights=const(B)
 assert j in (1,2,3) and h>F(j,5) and h<=F(4,5)
 smax=max(h,F(2,5))
 assert clo*expbounds((1-smax)*l)[0]>10**8
 # integer indices are rounded down in every log prime coordinate
 stop=int(h*(10*u));w=weights[:stop+1]+[0]*max(0,stop+1-len(weights))
 if stop>=120:
  w[120]+=up(S*(F(1,144)+F(1,162754)))
  for i in range(120,stop+1):w[i]+=ceil(S*(6*i*i-3*i+2),6*i**3)
 conv=w
 for _ in range(j-1):conv=[sum(conv[i]*w[t-i] for i in range(t+1)) for t in range(stop+1)]
 # R^2=k cN/(dL); includes all arithmetic-progression interval length remainders.
 def ab(k):
  D=rational_logupper(F(l)/(clo*k)); c=D+2*U
  if (1-smax)*l-D<=2*max(U,logupper(B)) or c>=(1-smax)*l: return None
  return 2*(1+F(k,l)+F(1,10**8)),1-c/l
 pairs=[b for k in (1,2,3,4,5) if (b:=ab(k)) is not None];assert pairs
 # Each component is monotone decreasing on [0,h] before truncation, so their minimum is also decreasing.
 def ker(t):
  x=F(t,10*u)
  return min(a*(h-x)/(b-x) for a,b in pairs)
 Ff=prod(range(6-j,6))*(h-F(j,5))
 total=F(sum(n*up(S*ker(i)) for i,n in enumerate(conv)),S**(j+1))/Ff
 sf=min(a/(b-F(2,5)) for a,b in pairs)
 square=sf*(F(sum(ceil(10**15,p*(p-1)) for p in primes if p>B),10**15)+F(1,n-1))+roundup(l/(clo*V*(expbounds(F(l,5))[0]-1)))
 return total,square
def run(l,u,B,h,j):
 a,s=upper(l,u,B,h,j);lo=cand(l,B)
 print(l,u,B,h,j,'candidate>',int(1000*lo),'bad<',up(1000*a),'square<',up(1000*s),'gap',round(float(lo-a-s),6),flush=True)
 return lo-a-s
if __name__=='__main__':
 for l,u,B,h,j in [(73,78,317,56,2),(78,83,317,57,2),(83,87,317,58,2),(87,89,317,58,2),(89,90,317,59,2),(90,91,317,59,2),(91,92,281,59,2)]:
  assert run(l,u,B,h,j)>0
