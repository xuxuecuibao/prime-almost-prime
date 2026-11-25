"""Four-factor estimates for the upper part of the finite interval, using the same weighted mean as fixed_platt_rz.py."""
from decimal import Decimal as D, ROUND_FLOOR
from math import factorial,prod
from functools import lru_cache
from fixed_platt_rz import c,box,H,K,one,primes,middle
def ev(x):return x.exp()
def floor(x):return int(x.to_integral_value(rounding=ROUND_FLOOR))
def par(l):
 if l<500: return 11,14,11,340
 if l<550: return 17,14,13,340
 if l<600: return 17,16,13,340
 if l<650: return 17,16,15,360
 if l<700: return 17,18,15,360
 if l<750: return 23,18,15,360
 if l<800: return 23,20,15,360
 if l<850: return 17,20,17,380
 if l<900: return 23,20,17,380
 if l<1000:return 23,22,17,380
 if l<1350:return 31,25,17,380
 if l<1550:return 31,31,17,380
 if l<1900:return 31,35,17,380
 if l<2150:return 31,39,17,380
 if l<2450:return 47,43,21,385
 if l<2650:return 53,45,21,390
 return 53,47,21,395
@lru_cache(None)
def geometry(b,vi):
 core=primes(b); v=D(vi)
 lam=v.ln()+D('.2615')+1/(2*v*v)-D('.5')-sum(one/p for p in core)
 lam+=sum(one/(p*(p-1)) for p in primes(10000) if p>b)+one/10000
 assert lam<3
 vt=ev(-lam/(1-one/(b-1)))
 de=vt*prod(1-one/(p-1) for p in core)
 return lam,vt,de,sum(D(p).ln() for p in core)
def row(l,u):
 b,vi,jb,mh=par(l)
 l,u,v,h=D(l),D(u),D(vi),D(mh)/1000
 lam,vt,de,q=geometry(b,vi); z=v.exp(); t=u**6; sc=D(27)
 kstar=floor(ev(sc)); dl=l/2-6*l.ln()
 assert q+jb*v<=dl and D(K)<kstar<t and z>K and ev(l/2)/l**6>t
 small=D('.006')/(l-3)
 eta=D(K)/c*(20*(D(K).ln()/2+2)*ev(-l/2)+((D(K)*(H+1)).ln()/2+2)*8/(H-1)*(2/(box*H))**24+2*ev(-l))
 outer=98*(4/u**D('2.5')+4/l**D('2.5')+18*l.sqrt()*ev(-l/12)+5*(1+u)*l**D('3.5')*ev(-l/6))/c
 def large(start):
  if start>=t:return outer
  return middle(l,u,start)/c+outer
 E=eta+large(K); Efar=large(kstar)
 U=max(one,lam-1,1-lam+lam**2/2)
 cap=one;n=1
 for p in primes(125):
  if n*p>kstar:break
  n*=p;cap*=one+one/(p-2)
 max_tail=floor(sc/D(b).ln())
 def mult(j):
  last=j-max_tail+1
  return min(U/de,cap+lam**last/factorial(last)/de) if last>=lam else U/de
 def gap(j):return lam**(j+1)/factorial(j+1)/vt
 s0=(1-small)*(1-gap(jb))-mult(jb)*E-U*Efar/de-D('1e-7')
 assert s0>0
 def ratio(j):
  return min(1+2*gap(j)/((1-(q+(j+1)*v)/l)*c*s0)+D('1e-7'),
             ((1+small)*(1+gap(j))+mult(j)*E+U*Efar/de+D('1e-7'))/s0)
 last=None;S=D(0)
 for j in range(20,1,-2):
  x=D(floor(1000*(D('.5')-(6*l.ln()+q+j*v)/l)))/1000
  if x<=v/l+D('.01'):continue
  x=min(x,h)
  assert q+j*v+x*l<=dl and q+(j+1)*v<=dl
  if last is None:
   f=h*((x*u/v).ln()+1/v**2+1/(z-1))-x+v/l
  else:
   f=h*(x/last).ln()-x+last+(h-last)*(1/(last*l)**2+1/(ev(last*l)-1))
  assert f>=0
  S+=f*ratio(j); last=x
  if x==h:break
 assert last==h
 eps=U*large(floor(z))/(de*s0)
 square=2/(D('.6')*(z-1)*c*de*s0)
 # Error in each comparison of adjacent sieve truncations via Brun--Titchmarsh.
 assert 4*u/(l**6*c*de*s0) <D('1e-7')
 # Prime powers, induced character replacements (including the p-principal part), large squares, omitted endpoints.
 assert (200*u**2*ev(-l/2)+4*u*ev(-l/5))/(c*de*s0)<D('1e-7')
 bad=(S+h*eps)/(5*h-1)+square+D('.000001')
 print(int(l),int(u),b,vi,jb,mh,f'{s0:.7f}',f'{bad:.7f}',flush=True)
 assert bad<D('.997')
 return bad
if __name__=='__main__':
 stops=[450,500,550,600,650,700,750,800,850,900,1000,1350,1550,1900,2150,2450,2650,2800]
 result=[]
 for a,b in zip(stops,stops[1:]):
  step=20 if a<1000 else 50
  for l in range(a,b,step):result.append(row(l,min(l+step,b)))
 print('max',max(result))
