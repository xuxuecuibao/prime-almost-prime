"""Four distinct factors on 2800 <= log N <= 84000. Python 3.11+.
Place next to upper_iv.py and its dependencies from the earlier certificate.
"""
from decimal import Decimal as D, ROUND_FLOOR, ROUND_CEILING
from math import factorial
from functools import lru_cache
from upper_iv import I,exp,log,mn,mx,primes,c,K,eta,middle,outer,one

def flo(x):return int(x.to_integral_value(rounding=ROUND_FLOOR))
def ceildiv(a,b):return (a+b-1)//b
@lru_cache(None)
def prime_consts(b):
 ps=primes(10000)
 q=I(0);s=I(0);cd=I(1);plus=1/I(10000)
 for p in ps:
  if p<=b:q+=log(p);s+=1/I(p);cd*=1-1/I(p-1)
  else:plus+=1/I(p*(p-1))
 return q,s,cd,plus
W=I(27);WI=flo(exp(W).lo)
assert flo(exp(W).hi)==WI
C=I(1);prod=1
for p in primes(125):
 if prod*p>WI:break
 prod*=p;C*=1+1/I(p-2)
def geo(b,w,v):
 q,sm,cd,plus=prime_consts(b)
 la=log(w)+I('.2615')+1/(2*w**2)-I('.5')-sm+plus
 lb=log(v/w)+1/(2*w**2)+1/(2*v**2)+1/(exp(w)-1)
 va=exp(-la/(1-1/I(b-1)));vb=exp(-lb/(1-1/(exp(w)-1)))
 assert la.lo>0 and lb.lo>0 and (w-log(b)).lo>0 and (v-w).lo>0
 return q,la,va,lb,vb,cd*va*vb
def gap(j,lam,mass):return lam**(j+1)/factorial(j+1)/mass
def omission(j,lam,mass,max_s):
 t=j-max_s+1
 if t>=lam.hi:return lam**t/factorial(t)/mass
 # Any individual intersection mass or omission of any alternating
 # partial sum versus its complete product is at most exp(lam).
 return exp(lam)/mass
def run(li,ui,b,wi,vi,ja,jb,mh):
 assert ja%2==0 and jb%2==0
 l,u,w,v,h=I(li),I(ui),I(wi),I(vi),I(mh)/100
 assert vi>=wi and (exp(v)-I(ui**6)).lo>0 # all detected conductors above T
 assert (exp(l/2)/l**6-I(ui**6)).lo>0
 q,la,va,lb,vb,de=geo(b,w,v)
 dl=l/2-6*log(l);base=q+(ja+1)*w+(jb+1)*v
 assert base.hi<dl.lo and (I('.5')-h-6/l).lo>0
 assert (5*h-1).lo>I('.5').hi and c.lo>D('.31')
 O=outer(l,u,True)
 E=eta(l)+middle(li,ui,K)/c+O
 Efar=middle(li,ui,WI)/c+O
 ga=gap(ja,la,va);gb=gap(jb,lb,vb)
 mas=flo((W/log(b)).hi);mbs=flo((W/w).hi)
 ma=omission(ja,la,va,mas);mb=omission(jb,lb,vb,mbs)
 Amax=C*((1+ma)*(1+omission(jb+1,lb,vb,mbs))+ma*(1+mb))
 un=3*exp(la+lb)
 s0=(1-I('.006')/(l-3))*(1-gb*(1+ga)-ga*(1+gb))-Amax*E-un*Efar/de-I('1e-7')
 assert s0.lo>D('.05') and (1-gb*(1+ga)-ga*(1+gb)).lo>0
 den=c*de*s0
 # For each level-D divisor or core-congruence term bound all exceptional
 # prime powers, small primes, nonreduced classes. Additional reserves cover
 # base/moment induction, deletion, and all p>N^.2 square divisors.
 assert (30*u/(l**6*den)).hi<D('1e-7')
 assert ((300*u**2*exp(-l/2)+4*u*exp(-l/5))/den).hi<D('1e-7')
 last=None;S=I(0)
 fixed=2*ga*exp(lb)/vb
 for n in range(jb,-1,-2):
  xi=flo((10000*(I('.5')-(6*log(l)+q+ja*w+n*v)/l)).lo)
  if I(xi).lo<=(10000*(v/l+I('.001'))).hi:continue
  t=mn(I(xi)/10000,h)
  if last is not None and t.lo<=last.hi:continue
  assert (q+ja*w+n*v+t*l).hi<dl.lo
  if last is None:
   F=h*(log(t*u/v)+1/v**2+1/(exp(v)-1))-t+v/l
  else:
   assert (last*l).lo>log(300).hi
   F=h*log(t/last)-t+last+(h-last)*(1/(last*l)**2+1/(exp(last*l)-1))
  assert F.lo>0
  gn=gap(n,lb,vb)
  denA=1-(q+(ja+1)*w+n*v)/l
  denB=1-((n+1)*v+2*(w+2)+log(l/c))/l
  # Montgomery with R²=cN/(dL): R>e^w and log R>U_B.
  assert denA.lo>0 and denB.lo>0
  assert (c*exp(l-(n+1)*v)).lo>10**8
  assert (2*(1+1/l+I('1e-8'))).hi<D('2.003')
  rb=1+(fixed/(c*denA)+I('2.003')*gn/denB)/s0+I('1e-7')
  mbn=omission(n,lb,vb,mbs)
  rc=((1+I('.006')/(l-3))*(1+ga)*(1+gn)+C*(1+ma)*(1+mbn)*E+exp(la+lb)*Efar/de+I('1e-7'))/s0
  S+=F*mn(rb,rc);last=t
  if last.lo==h.lo and last.hi==h.hi:break
 assert last is not None and last.lo==h.lo and last.hi==h.hi
 eps=exp(la+lb)*O/(de*s0)
 square=2/(I('.6')*(exp(v)-1)*den)
 bad=(S+h*eps)/(5*h-1)+square+I('1e-6')
 assert bad.hi<D('1'),(li,ui,bad.hi)
 return s0.lo,bad.hi
def rows():
 # Boundaries simply choose sieve parameters; rows are contiguous.
 regimes=[(2800,4000,53,38,5,10),(4000,6500,53,38,5,12),
          (6500,13000,101,38,5,14),(13000,25000,101,42,5,16),
          (25000,55000,101,50,5,20),(55000,84000,503,60,5,20)]
 for lo,hi,b,ir,rat,ja in regimes:
  l=lo
  while l<hi:
   u=min(hi,l*104//100)
   v=ceildiv(l,ir);w=ceildiv(v,rat)
   # floor uses only an approximation to choose an even degree.
   # run() rigorously checks every actual logarithmic support.
   from math import log as fl
   jb=2*int((l/2-6*fl(l)-b-(ja+1)*w-v)//(2*v))
   yield (l,u,b,w,v,ja,jb,37)
   l=u
if __name__=='__main__':
 last=2800;worst=D(0)
 for row in rows():
  assert row[0]==last;last=row[1]
  s0,bd=run(*row);worst=max(worst,bd)
  print(*row,'s0 >=',s0.quantize(D('.00000001'),rounding=ROUND_FLOOR),'bad <=',bd.quantize(D('.00000001'),rounding=ROUND_CEILING),flush=True)
 assert last==84000;print('worst',worst)
