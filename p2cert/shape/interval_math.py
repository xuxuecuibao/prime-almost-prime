"""55-digit outward-rounded interval operations from the attached p4_fd_verify.py.
Only the arithmetic layer, odd prime enumeration and weight constants are required.
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
