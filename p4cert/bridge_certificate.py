"""Interval operations and the grouped Type II estimate imported by upper_iv.py.
"""
from decimal import Decimal, getcontext, ROUND_CEILING
from math import isqrt
from pathlib import Path
getcontext().prec = 55
D=Decimal

class I:
    def __init__(self, x, y=None):
        self.lo = D(str(x))
        self.hi = self.lo if y is None else D(str(y))
    def __add__(self,b):
        b=ival(b)
        return I((self.lo+b.lo).next_minus(),(self.hi+b.hi).next_plus())
    __radd__=__add__
    def __neg__(self): return I(self.hi.copy_negate(),self.lo.copy_negate())
    def __sub__(self,b):return self+-ival(b)
    def __rsub__(self,b):return ival(b)+-self
    def __mul__(self,b):
        b=ival(b); v=[x*y for x in (self.lo,self.hi) for y in (b.lo,b.hi)]
        return I(min(v).next_minus(),max(v).next_plus())
    __rmul__=__mul__
    def __truediv__(self,b):
        b=ival(b);assert not b.lo<=0<=b.hi
        return self*I((1/b.hi).next_minus(),(1/b.lo).next_plus())
    def __rtruediv__(self,b):return ival(b)/self
    def __pow__(self,n):
        assert isinstance(n,int) and n>=0 and self.lo>=0
        r=I(1);x=self
        while n:
            if n%2:r=r*x
            n//=2;x=x*x
        return r

def ival(x): return x if isinstance(x,I) else I(x)
def exp(x):
    x=ival(x); return I(x.lo.exp().next_minus(),x.hi.exp().next_plus())
def log(x):
    x=ival(x); return I(x.lo.ln().next_minus(),x.hi.ln().next_plus())
def sqrt(x):
    x=ival(x); return I(x.lo.sqrt().next_minus(),x.hi.sqrt().next_plus())
def mn(x,y):
    x,y=ival(x),ival(y);return I(min(x.lo,y.lo),min(x.hi,y.hi))
def mx(x,y):
    x,y=ival(x),ival(y);return I(max(x.lo,y.lo),max(x.hi,y.hi))
def ceiling_upper(x):return int(x.hi.to_integral_value(rounding=ROUND_CEILING))
def primes(n):
    p=bytearray([1])*(n+1)
    for i in range(2,isqrt(n)+1):
        if p[i]: p[i*i:n+1:i]=bytes((n-i*i)//i+1)
    return [i for i in range(3,n+1,2) if p[i]]

K=400000; H=I('293.75'); d=I('5.4')/H; c=(1-exp(-1))**2*((1-exp(-d))/d)**24
ROWS=[(159,170,103,7,2,252),(170,185,151,7,2,264),
      (185,200,211,11,2,272),(200,212,283,13,2,278)]

def branch(l,u,q,choice):
    a,t,s,nu=map(I,choice)
    if a.lo==0:
        return 98*(u**3)*sqrt(u)*(4+2*exp(2*q-l/2)+6*exp(3*q/2-l/3)+5*exp(q-l/6))
    v=mx(a*q,log(101)); v0=t*v; lr=log(nu); lu=l-s*q-v0-2*log(u)
    assert v.lo>=log(100).hi and lu.lo>=1
    assert (s*q+v0+2*log(u)).lo>max(v.hi,(2+24*d).hi)
    ty=(I('2.08')*exp((3-s)*q)/u**2+(2*u+4)*exp(v0+3*q-l))*mn(1,(1+q)*exp(-q/2))
    Y=lr+s*q+2*log(u)+v0; g=nu/(nu-1)
    J=1+ceiling_upper((Y-v)/lr)
    z=lu-lr
    X=(nu-1+I('.005')*(nu+1)/z+2*exp(-z/2)*(z+1) if z.lo>30 else I('1.04')*nu)
    G=exp(1+I('.578')/Y)*Y/((t-1)**2*v)*(I('1.084')*(t+1)+
          (I('1.301')-I('.318')*exp(-v/Y)+(I('1.301')-I('.318')*exp(-v0/Y))*t**2)*v/Y)
    grouped=(nu-1)*J+g*nu*(exp(2*q-lu)+exp(-lu))+g*(nu-1)*exp(2*q-v)+J*(exp(4*q-l)+exp(2*q-l))
    return ty+sqrt(X*u*grouped*g*G)

def supports(l,u,logD0,logD1,eta):
    recipes=Path(__file__).with_name('bridge_choices.csv').read_text().splitlines()
    idx=0; acc=I(0); zprev=1/I(K); E=eta if logD0.hi<=log(K).lo else None
    for row in recipes:
        if row.startswith('#'):continue
        a,b,start,stop,*choice=row.split(',')
        if int(a)!=l.lo or int(b)!=u.hi:continue
        assert int(start)==idx+1
        for idx in range(int(start),int(stop)+1):
            q=log(K)+idx*I('.12'); z=exp(-q); A=branch(l,u,q,choice)
            acc+=A*(zprev-z); zprev=z
            if E is None and q.lo>=logD0.hi:E=eta+(acc+A*z)/c
            if q.lo>=logD1.hi:
                assert E is not None
                return E,eta+(acc+A*z)/c
    raise ArithmeticError('incomplete support')

def run(l0,u0,B,b,j,mh):
    l=I(l0);u=I(u0);h=I(mh)/1000;y=log(B);allp=primes(10000)
    qc=sum((log(p) for p in allp if p<=b),I(0))
    v=I(1);delta=I(1);T=[I(1)]+[I(0)]*(j+2)
    for p in allp:
        if p>B:break
        a=1/I(p-1);delta*=1-a
        if p<=b:continue
        v*=1-a
        for k in range(j+2,0,-1):T[k]+=a*T[k-1]
    assert T[1].hi<2
    sq=sum((1/I(p*(p-1)) for p in allp if p>B),1/I(10000))
    rho=I(0);t=I(0)
    for p in allp:
        if B<p<=300:
            t+=1/I(p-1);rho=mx(rho,t-log(log(p)/y))
    f=log(300);rho=mx(rho,t-log(f/y)+1/f**2+1/I(299))
    logD0=qc+(j+1)*y;logD1=qc+j*y+h*u
    assert logD0.hi<l.lo and (h*l).lo>y.hi
    eta=K/c*(20*(log(K)/2+2)*exp(-l/2)+(log(K*(H+1))/2+2)*8/(H-1)*(2/(d*H))**24+2*exp(-l))
    E,P=supports(l,u,logD0,logD1,eta)
    m1=(1-exp(-2))/c*I('.0024')/(l-2-24*d)
    safety=I('1e-7')
    s0=(1-m1)*(1-T[j+2]/v)-E/delta-safety
    assert s0.lo>0
    R=mn(((1+m1)*(1+T[j+1]/v)+E/delta+safety)/s0,
         1+2*T[j+1]/(c*(1-logD0/l)*s0*v)+safety)
    W=h*(log(h*u/y)-1+rho)+y/l
    bound=(R*W+h*P/(s0*delta))/(6*h-1)+2*sq/(I('.6')*s0*delta*c)+I('1e-6')
    # Reserve checks, normalized conservatively by the lower bound on S.
    mass=c*delta*s0
    induction=(b+j+2)*u*(exp(logD0-l)+2*exp(logD1-l))/mass
    primepowBT=I(B)**(j+1)*(u*exp(-l/2)+(b+1)*u*exp(-l))/mass
    large_square=u*(1/(exp(l/5)-1)+exp(-l/2))
    deletions=u*(exp(-l/2)+(B+3+u)*exp(-l))
    assert induction.hi<D('1e-12') and primepowBT.hi<D('1e-12')
    assert ((large_square+deletions)/mass).hi<D('1e-7')
    assert bound.hi<D('.99')
    return s0.lo,E.hi,P.hi,bound.hi

if __name__=='__main__':
    for r in ROWS:
        print(*r,*(format(x,'.10f') for x in run(*r)),flush=True)
