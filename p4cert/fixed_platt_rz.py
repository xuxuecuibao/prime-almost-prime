"""Finite certificates for a prime plus a P5. Run with Python 3.11+.
Parameter choice uses math (affecting only speed/strength). All selected
bounds and comparisons are recalculated at 45 decimal digits.
"""
from decimal import Decimal as D, getcontext, ROUND_CEILING, ROUND_FLOOR
from functools import lru_cache
from math import ceil, exp, log, sqrt, prod, isqrt, factorial
getcontext().prec=45
K=400000
H=D('293.75')
box=D('5.4')/H
rbox=24
one=D(1)
c=(1-(-one).exp())**2*((1-(-box).exp())/box)**rbox

def branches(l,u,q,choose=None):
    highprec=choose is not None
    if highprec:
        expf=lambda x:x.exp(); logf=lambda x:x.ln(); sq=lambda x:x.sqrt()
        ceilf=lambda x:int(x.to_integral_value(rounding=ROUND_CEILING)); f=lambda s:D(str(s))
    else:
        expf=exp;logf=log;sq=sqrt;ceilf=ceil;f=float
    ah=98*u**f(3.5)*(4+2*expf(2*q-f(.5)*l)+6*expf(f(1.5)*q-l/3)+5*expf(q-l/6))
    if choose==():return ah
    best,params=ah,()
    for d in ([choose[0]] if highprec else (1.25,1.5,1.75,2,2.25,2.5)):
        v=f(d)*q
        if v<logf(f(100)):continue
        for ta in ([choose[1]] if highprec else (1.25,1.5,1.75,2,2.5,3,4,6)):
            t=f(ta);v0=t*v
            for rr in ([choose[2]] if highprec else (1.5,2,2.25,2.5,2.75,3)):
                r=f(rr);lu=l-r*q-v0-2*logf(u)
                if lu<1:continue
                ty=(f(2.08)*expf((3-r)*q)/u**2+(2*u+4)*expf(v0+3*q-l))*min(f(1),(1+q)*expf(-q/2))
                for n in ([choose[3]] if highprec else (1.08,1.2,1.5,2)):
                    nu=f(n);lr=logf(nu);Y=lr+r*q+2*logf(u)+v0
                    J=1+ceilf((Y-v)/lr)
                    x=(nu-1+f(.005)*(nu+1)/(lu-lr)+2*expf(-(lu-lr)/2)*(lu-lr+1) if lu-lr>30 else f(1.04)*nu)
                    mf=nu-1+nu*expf(2*q-lu)+nu*expf(-lu)
                    E=expf(1+f(.578)/Y)*(f(1.084)*(t+1)+(f(1.301)-f(.318)*expf(-v/Y)+(f(1.301)-f(.318)*expf(-v0/Y))*t**2)*v/Y)*Y/((t-1)**2*v)
                    dy=1/(1-1/nu)
                    bound=ty+sq(x*mf*u*(J+dy*expf(2*q-v))*dy*E)
                    if highprec:return bound
                    if bound<best:best,params=bound,(d,ta,rr,n)
    if highprec:raise ArithmeticError('invalid selected parameters')
    return params

@lru_cache(None)
def middle(l,u,start):
    l,u,start=map(D,(l,u,start));T=u**6; y=start; total=D(0)
    while y<T:
        z=min(T,y*D('1.25')); q=z.ln()
        choice=branches(float(l),float(u),float(q))
        bound=branches(l,u,q,choice)
        total+=bound*(1/y-1/z);y=z
    return total+bound/T

@lru_cache(None)
def primes(N):
    a=bytearray(b'\x01')*(N+1);a[:2]=b'\0\0'
    for p in range(2,isqrt(N)+1):
        if a[p]:a[p*p::p]=b'\0'*((N-p*p)//p+1)
    return tuple(p for p in range(3,N+1,2) if a[p])

@lru_cache(None)
def sieve_numbers(B,b,j):
    P=min(B,10000);y=D(B).ln(); vv=D(P);v0=vv.ln();ps=primes(P);tail=[p for p in ps if p>b]
    v=prod(1-one/(p-1) for p in tail);delta=v*prod(1-one/(p-1) for p in ps if p<=b)
    X=[one]+[D(0)]*(j+2)
    for p in tail:
        a=one/(p-1)
        for k in range(j+2,0,-1): X[k]+=a*X[k-1]
    if B>P:
        A=(y/v0).ln()+1/(2*y*y)+1/(2*v0*v0)+1/(vv-1)
        tailprod=(-A/(1-one/(vv-1))).exp();v*=tailprod;delta*=tailprod
        X=[sum(X[n-k]*(A**k if k else 1)/factorial(k) for k in range(n+1)) for n in range(j+3)]
    sq=(sum(one/(p*(p-1)) for p in primes(10000) if p>B)+one/10000 if B<10000 else one/(B-1))
    if B<300:
        t=D(0);rho=D(0)
        for p in primes(300):
            if p>B:
                t+=one/(p-1);rho=max(rho,t-(D(p).ln()/y).ln())
        f=D(300).ln();rho=max(rho,t-(f/y).ln()+1/(f*f)+one/299)
    else:rho=1/y**2+one/(B-1)
    q=sum(D(p).ln() for p in ps if p<=b)
    return q,v,delta,X,sq,rho

@lru_cache(None)
def arithmetic(l,u,B,b,j,mh=0):
    l,u=map(D,(l,u));q,v,d,X,sq,rho=sieve_numbers(B,b,j)
    y=D(B).ln(); h=(D(mh)/1000 if mh else ((D('.5')-(q+j*y+6*l.ln())/l)*1000).to_integral_value(rounding=ROUND_FLOOR)/1000)
    assert X[1]<2 and 6*h>1 and q+j*y+h*l<=l/2-6*l.ln() and q+(j+1)*y<=l/2-6*l.ln()
    m1=(1-(-D(2)).exp())/c*D('.0024')/(l-2-rbox*box)
    eta=K*(20*(D(K).ln()/2+2)*(-l/2).exp()+((D(K)*(H+1)).ln()/2+2)*8/(H-1)*(2/(box*H))**rbox+2*(-l).exp())/c
    outer=98*(4/u**D('2.5')+4/l**D('2.5')+18*l.sqrt()*(-l/12).exp()+5*(1+u/2-6*u.ln()-6*l.ln())*l**D('3.5')*(-l/6).exp())/c
    E=eta+middle(l,u,K)/c+outer
    Perror=(eta if B<K else 0)+middle(l,u,max(B,K))/c+outer
    safety=D('1e-7');s0=(1-m1)*(1-X[j+2]/v)-E/d-safety
    assert s0>0
    R=min(((1+m1)*(1+X[j+1]/v)+E/d+safety)/s0,1+2*X[j+1]/(1-(q+(j+1)*y)/l)/c/s0/v+safety)
    W=h*((h*u/y).ln()-1+rho)+y/l
    bad=(R*W+h*Perror/(s0*d))/(6*h-1)+2*sq/(D('.6')*s0*d*c)+D('1e-6')
    denom=c*d*s0
    induction=4*(b+j+2)*(u+1)*(-l/2).exp()/denom
    diff_pp=D(B)**(j+1)*u*(-l/2).exp()/denom
    deletions=((D(2)*(-l/2).exp()+(B+3+u)*(-l).exp())*u+u/((l/5).exp()-1))/denom
    assert induction<D('1e-12') and diff_pp<D('1e-9') and deletions<D('1e-7')
    return h,s0,E,Perror,bad

ROWS=[(212,217,110,5,4),(217,221,110,5,4),(221,230,170,5,4),(230,242,170,5,4),
      (242,266,286,7,4),(266,300,286,7,4),(300,340,500,11,4),(340,368,1000,13,4),
      (368,460,700,13,4),(460,570,5000,7,6),(570,680,20000,13,6),(680,820,50000,17,6),
      (820,1000,100000,23,6),(1000,1200,300000,7,8),(1200,1500,1000000,13,8),
      (1500,1800,1000000,13,8),(1800,2100,10000000,23,8),(2100,2500,100000000,31,8),
      (2500,2800,100000000,31,8)]
if __name__=='__main__':
    for args in ROWS:
        h,s0,E,Perr,bad=arithmetic(*args)
        print(*args,*(f'{w:.7f}' for w in (h,s0,E,Perr,bad)),flush=True)
        assert bad<D('.997')
