"""Directed enclosures on every upper row (212 through 2800).
Place beside bridge_certificate.py, four_variable.py and four_bands.py.
Any double-precision work selects a candidate; no floating value is used as
an upper or lower bound.  Uses Python's standard library.
"""
from math import exp as fe, log as fl, sqrt as fs, ceil, factorial
from decimal import Decimal as D, ROUND_FLOOR, getcontext
from functools import lru_cache
from four_variable import ROWS
from four_bands import par
from bridge_certificate import I,ival,exp,log,sqrt,mn,mx,primes,branch,K,H,d,c
getcontext().prec=55
one=I(1)
def floor(x): return int(x.to_integral_value(rounding=ROUND_FLOOR))
def exactfloor(x):
    assert floor(x.lo)==floor(x.hi)
    return floor(x.lo)
def cutoff(x,k): return I(floor((x*k).lo))/k
def point(x): return x.lo if x.lo==x.hi else (x.hi+x.lo)/2

def choose(l,u,q):
    # Small menu only. branch() reruns the winning parameter choice with
    # outward rounding, using its independently proved grouped bound.
    best=98*u**3.5*(4+2*fe(2*q-l/2)+6*fe(1.5*q-l/3)+5*fe(q-l/6))
    choice=(0,0,0,0)
    for dd in (1.75,2.):
      for t in (1.25,1.5,1.75,2.,2.5,3.,4.):
       for r in (1.5,2.,2.25):
        v=dd*q; v0=t*v; lu=l-r*q-v0-2*fl(u)
        if lu<1.01: continue
        ty=(2.08*fe((3-r)*q)/u**2+(2*u+4)*fe(v0+3*q-l))*min(1.,(1+q)*fe(-q/2))
        if ty>best: continue
        for n in ((1.08,1.2,1.5) if t<=1.5 else (1.08,)):
          lr=fl(n);Y=lr+r*q+2*fl(u)+v0;g=n/(n-1);J=1+ceil((Y-v)/lr);z=lu-lr
          xx=n-1+.005*(n+1)/z+2*fe(-z/2)*(z+1) if z>30 else 1.04*n
          gg=fe(1+.578/Y)*Y/((t-1)**2*v)*(1.084*(t+1)+(1.301-.318*fe(-v/Y)+(1.301-.318*fe(-v0/Y))*t**2)*v/Y)
          group=(n-1)*J+g*n*(fe(2*q-lu)+fe(-lu))+g*(n-1)*fe(2*q-v)+J*(fe(4*q-l)+fe(2*q-l))
          z=ty+fs(xx*u*group*g*gg)
          if z<best: best=z; choice=(dd,t,r,n)
    return choice
@lru_cache(None)
def middle(li,ui,start):
    l,u,T,y=I(li),I(ui),I(ui**6),I(start); acc=I(0)
    assert start<ui**6
    while y.hi<T.lo:
       z=mn(T,y*I('1.25'));q=log(z)
       ch=choose(li,ui,float(point(q)))
       A=branch(l,u,q,ch)
       acc+=A*(1/y-1/z);y=z
    assert y.lo==T.lo and y.hi==T.hi
    return acc+A/T
def eta(l):
    return K/c*(20*(log(K)/2+2)*exp(-l/2)+(log(K*(H+1))/2+2)*8/(H-1)*(2/(d*H))**24+2*exp(-l))
def outer(l,u,coarse):
    coef=1+u if coarse else 1+u/2-6*log(u)-6*log(l)
    assert coef.lo>0
    return 98*(4/(u**2*sqrt(u))+4/(l**2*sqrt(l))+18*sqrt(l)*exp(-l/12)+5*coef*l**3*sqrt(l)*exp(-l/6))/c

@lru_cache(None)
def sieve(B,b,j):
    P=min(B,10000);y=log(B);vp=log(P);v=one;delta=one
    X=[one]+[I(0) for _ in range(j+2)]
    ps=primes(10000)
    for p in ps:
        if p>P: break
        a=1/I(p-1);delta*=1-a
        if p<=b:continue
        v*=1-a
        for k in range(j+2,0,-1):X[k]+=a*X[k-1]
    if B>P:
        A=log(y/vp)+1/(2*y**2)+1/(2*vp**2)+1/I(P-1)
        w=exp(-A/(1-1/I(P-1)));v*=w;delta*=w
        X=[sum((X[n-k]*A**k/factorial(k) for k in range(n+1)),I(0)) for n in range(j+3)]
    sq=sum((1/I(p*(p-1)) for p in ps if p>B),1/I(10000)) if B<10000 else 1/I(B-1)
    assert B>=300
    rho=1/y**2+1/I(B-1)
    q=sum((log(p) for p in ps if p<=b),I(0))
    return q,v,delta,X,sq,rho
def variable(args):
    li,ui,B,b,jb,mh=args
    l,u=I(li),I(ui);h=I(mh)/1000;y=log(B);q,v,de,X,sq,rho=sieve(B,b,jb)
    assert (5*h).lo>1 and X[1].hi<2
    mi=(1-exp(-2))/c*I('.0024')/(l-2-24*d)
    E=eta(l)+middle(li,ui,K)/c+outer(l,u,False)
    P=E  # all these B are below K
    assert B<K
    s0=(1-mi)*(1-X[jb+2]/v)-E/de-I('1e-7')
    level=l/2-6*log(l);assert s0.lo>0 and (q+(jb+1)*y).hi<level.lo
    den=c*s0*de
    assert (I(B)**(jb+1)*u*exp(-l/2)/den).hi<D('1e-9')
    assert ((4*(b+jb+2)*(u+2)*exp(-l/2)+2*u/(exp(l/5)-1))/den).hi<D('1e-8')
    S=I(0);last=None
    for j in range(jb,-1,-2):
       t=mn(h,cutoff((level-q-j*y)/l,10**8))
       if t.lo<=(y/l).hi or (last is not None and t.lo<=last.hi):continue
       assert (q+j*y+t*l).hi<=level.lo and (I('.5')-t-6/l).lo>0
       if last is None:F=h*(log(t*u/y)+rho)-t+y/l
       else:
           a=last*l;F=h*log(t/last)-t+last+(h-last)*(1/a**2+1/(exp(a)-1))
       assert F.lo>=0
       g=X[j+1]/v
       R=mn(((1+mi)*(1+g)+E/de+I('1e-7'))/s0,1+2*g/((1-(q+(j+1)*y)/l)*c*s0)+I('1e-7'))
       S+=R*F;last=t
       if t.lo==h.lo and t.hi==h.hi:break
    assert last.lo==h.lo and last.hi==h.hi
    bad=(S+h*P/(s0*de))/(5*h-1)+2*sq/(I('.6')*den)+I('1e-6')
    assert bad.hi<D('.99')
    return s0.lo,bad.hi

@lru_cache(None)
def geometry(b,vi):
    core=primes(b);v=I(vi);ps=primes(10000)
    lam=log(v)+I('.2615')+1/(2*v**2)-I('.5')-sum((1/I(p) for p in core),I(0))
    lam+=sum((1/I(p*(p-1)) for p in ps if p>b),1/I(10000))
    assert lam.hi<3
    vt=exp(-lam/(1-1/I(b-1)));de=vt
    for p in core:de*=1-1/I(p-1)
    return lam,vt,de,sum((log(p) for p in core),I(0))

def bands(li,ui):
    b,vi,jb,mh=par(li); l,u,v,h=I(li),I(ui),I(vi),I(mh)/1000
    lam,vt,de,q=geometry(b,vi);z=exp(v);kstar=exactfloor(exp(27));dl=l/2-6*log(l)
    assert (q+jb*v).hi<dl.lo and K<kstar<ui**6 and z.lo>K and (exp(l/2)/l**6).lo>ui**6
    small=I('.006')/(l-3);out=outer(l,u,True)
    def large(start):return out if start>=ui**6 else middle(li,ui,start)/c+out
    E=eta(l)+large(K);Efar=large(kstar)
    U=mx(one,mx(lam-1,1-lam+lam**2/2))
    cap=one;n=1
    for p in primes(125):
        if n*p>kstar:break
        n*=p; cap*=1+1/I(p-2)
    mt=exactfloor(27/log(b))
    def mult(j):
       r=j-mt+1
       return mn(U/de,cap+lam**r/factorial(r)/de) if r>=lam.hi else U/de
    def gap(j):return lam**(j+1)/factorial(j+1)/vt
    s0=(1-small)*(1-gap(jb))-mult(jb)*E-U*Efar/de-I('1e-7')
    assert s0.lo>0 and ((I('.5')-h)*l).lo>6
    def ratio(j):
       return mn(1+2*gap(j)/((1-(q+(j+1)*v)/l)*c*s0)+I('1e-7'),((1+small)*(1+gap(j))+mult(j)*E+U*Efar/de+I('1e-7'))/s0)
    last=None;S=I(0)
    for j in range(20,1,-2):
       x=cutoff(I('.5')-(6*log(l)+q+j*v)/l,1000)
       if x.lo<=(v/l+I('.01')).hi:continue
       x=mn(h,x)
       assert x.lo>(v/l if last is None else last).hi
       assert (q+j*v+x*l).hi<dl.lo and (q+(j+1)*v).hi<dl.lo
       if last is None:F=h*(log(x*u/v)+1/v**2+1/(z-1))-x+v/l
       else:F=h*log(x/last)-x+last+(h-last)*(1/(last*l)**2+1/(exp(last*l)-1))
       assert F.lo>=0
       S+=F*ratio(j);last=x
       if x.lo==h.lo and x.hi==h.hi:break
    assert last.lo==h.lo and last.hi==h.hi
    eps=U*large(exactfloor(z))/(de*s0)
    square=2/(I('.6')*(z-1)*c*de*s0)
    assert (4*u/(l**6*c*de*s0)).hi<D('1e-7')
    assert ((200*u**2*exp(-l/2)+4*u*exp(-l/5))/(c*de*s0)).hi<D('1e-7')
    bad=(S+h*eps)/(5*h-1)+square+I('.000001')
    assert bad.hi<D('.997')
    return s0.lo,bad.hi
def rows():
    st=[450,500,550,600,650,700,750,800,850,900,1000,1350,1550,1900,2150,2450,2650,2800]
    return [(l,min(l+(20 if a<1000 else 50),b)) for a,b in zip(st,st[1:]) for l in range(a,b,20 if a<1000 else 50)]
if __name__=='__main__':
    from sys import argv
    last=212;results=[]
    def report(li,ui,res):
        global last
        assert li<=last<ui; last=ui;results.append(res[1])
        print(li,ui,'s0 >=',res[0].quantize(D('.00000001'),rounding=ROUND_FLOOR),'bad <=',res[1].quantize(D('.00000001'),rounding='ROUND_CEILING'),flush=True)
    if 'bands' not in argv:
       for a in ROWS:report(a[0],a[1],variable(a))
       print('first maximum',max(results),flush=True)
    if 'variable' not in argv:
       last=450;results=[]
       for li,ui in rows():report(li,ui,bands(li,ui))
       assert last==2800
       print('second maximum',max(results),flush=True)
