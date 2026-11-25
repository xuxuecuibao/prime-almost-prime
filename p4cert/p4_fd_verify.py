"""Outward-rounded certificate for variable-degree four-factor subintervals (standard library).
The parameter rows choose admissible upper bounds only.
"""
from decimal import Decimal, getcontext, ROUND_CEILING
from math import isqrt
from parameter_tables import FINITE_DIFFERENCE_CHOICES, iter_parameter_rows
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
from functools import lru_cache
def dcon(n):
    z=[1,2]+[4*250**(k-2) for k in range(2,n+1)]
    row=[1]
    for t in range(n):
        row=[(t*row[k] if k<len(row) else 0)+(row[k-1] if k else 0) for k in range(len(row)+1)]
    return I('.501')*4/I(4**n)*exp(I('2.45')*(n-1))*sum(x*y for x,y in zip(z,row))
C=[(i,dcon(i)) for i in (2,3,4,5,6)]
def branch(l,u,q,choice):
    if not choice:
        return 98*u**3*sqrt(u)*(4+2*exp(2*q-l/2)+6*exp(3*q/2-l/3)+5*exp(q-l/6))
    a,t,s,nu=map(I,choice)
    v=mx(a*q,log(101)); v0=t*v; lr=log(nu); lu=l-s*q-v0-2*log(u)
    assert v.lo>=log(100).hi and lu.lo>=1
    assert (s*q+v0+2*log(u)).lo>max(v.hi,(2+24*d).hi)
    ty=(I('2.08')*exp((3-s)*q)/u**2+(2*u+4)*exp(v0+3*q-l))*mn(1,(1+q)*exp(-q/2))
    for k,cn in C:
        exponent=q*(I(k)+I('1.5'))
        d1=log(I('1.04'))+exponent-k*s*q-2*k*log(u)
        d2=log(u+k)+k*(v0-l)+exponent
        ty=mn(ty,cn*(exp(d1)+exp(d2)))
    Y=lr+s*q+2*log(u)+v0; g=nu/(nu-1); J=1+ceiling_upper((Y-v)/lr)
    z=lu-lr
    X=nu-1+I('.005')*(nu+1)/z+2*exp(-z/2)*(z+1) if z.lo>30 else I('1.04')*nu
    G=exp(1+I('.578')/Y)*Y/((t-1)**2*v)*(I('1.084')*(t+1)+(I('1.301')-I('.318')*exp(-v/Y)+(I('1.301')-I('.318')*exp(-v0/Y))*t**2)*v/Y)
    grouped=(nu-1)*J+g*nu*(exp(2*q-lu)+exp(-lu))+g*(nu-1)*exp(2*q-v)+J*(exp(4*q-l)+exp(2*q-l))
    return ty+sqrt(X*u*grouped*g*G)
alpha=I('4.1');GSTEP=I('1.25')
def cum(l,u,idx,lines):
    top=u/2-alpha*log(u);start=log(K);last=1/I(K);acc=I(0)
    for i,ch in enumerate(lines,1):
        q=mn(start+i*log(GSTEP),top)
        A=branch(l,u,q,ch)
        x=exp(-q);acc+=A*(last-x);last=x
        if (start+i*log(GSTEP)).lo>=top.hi:return (acc+A*x)/c
    raise AssertionError('incomplete conductors')
P=primes(20000)
@lru_cache(None)
def lprime(p):return log(p)
@lru_cache(None)
def local(B,b,j):
    y=log(B);v=I(1);delta=I(1);T=[I(1)]+[I(0)]*(j+2)
    for p in P:
        if p>B:break
        a=1/I(p-1);delta*=1-a
        if p<=b:continue
        v*=1-a
        for k in range(j+2,0,-1):T[k]+=a*T[k-1]
    assert T[1].hi<2
    sq=sum((1/I(p*(p-1)) for p in P if p>B),1/I(20000))
    qc=sum((lprime(p) for p in P if p<=b),I(0))
    # Rigorous cumulative prime mass through 10,000 followed by RS bounds.
    rho=I(0);t=I(0)
    for p in P:
        if B<p<=10000:
            t+=1/I(p-1);rho=mx(rho,t-log(lprime(p)/y))
    f=log(10000);rho=mx(rho,t-log(f/y)+1/f**2+1/I(9999))
    return y,v,delta,T,sq,qc,rho
def run(row,lines):
    L,U,B,b,j,mh=row;l=I(L);u=I(U);h=I(mh)/100
    y,v,delta,T,sq,qc,rho=local(B,b,j)
    assert (qc+(j+1)*y).hi< (l/2-alpha*log(l)).lo
    eta=K/c*(20*(log(K)/2+2)*exp(-l/2)+(log(K*(H+1))/2+2)*8/(H-1)*(2/(d*H))**24+2*exp(-l))
    E=eta+cum(l,u,0,lines)
    m1=(1-exp(-2))/c*I('.0024')/(l-2-24*d)
    s0=(1-m1)*(1-T[j+2]/v)-E/delta-I('1e-7')
    assert s0.lo>0
    den=c*delta*s0
    ind=4*(b+j+2)*(u+1)*exp(-l/2)/den
    assert ind.hi<D('1e-11')
    last=I(0);wm=I(0)
    for n in range(j,-1,-2):
        cutoff=I('.5')-(qc+n*y+alpha*log(l))/l
        if cutoff.hi<=(y/l).lo:continue
        assert cutoff.lo>(y/l).hi
        if cutoff.lo>=h.hi:tau=h;complete=True
        else:
            assert cutoff.hi<h.lo
            tau=cutoff;complete=False
        if last.hi==0:
            w=tau*(log(tau*u/y)-1+rho)+y/l+(h-tau)*(log(tau*u/y)+rho)
        else:
            Y=last*l
            assert Y.lo>log(300).hi
            w=h*log(tau/last)-(tau-last)+(h-last)*(1/Y**2+1/(exp(Y)-1))
        pp=I(B)**(n+1)*(u*exp(-l/2)+(b+1)*u*exp(-l))/den
        r1=((1+m1)*(1+T[n+1]/v)+E/delta+I('1e-7'))/s0
        r2=1+2*T[n+1]/(c*(1-(qc+(n+1)*y)/l)*s0*v)+pp
        wm+=mn(r1,r2)*w;last=tau
        if complete:break
    assert complete and (5*h-1).lo>0 and h.hi<I('.5').lo
    delete=u*(2*exp(-l/2)+(B+3+u)*exp(-l)+1/(exp(l/5)-1))/den
    assert delete.hi<D('0.00001')
    bad=(wm+h*(E+ind)/(s0*delta))/(5*h-1)+2*sq/(I('.6')*den)+delete+I('.000001')
    print(L,U,B,b,j,mh,'error',format(E.hi,'.8f'),'bad',format(bad.hi,'.8f'),flush=True)
    assert bad.hi<D('.977714')
    return bad.hi
def main():
    for row, choices in iter_parameter_rows(FINITE_DIFFERENCE_CHOICES):
        run(row, choices)

if __name__=='__main__':main()
