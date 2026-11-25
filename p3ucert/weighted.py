"""Uniform numerical inequalities in the P3 proof, for log N >= 23500.
Run python3 -B weighted.py. Only the standard library and interval_math.py are used.
"""
from interval_math import I, log, exp, sqrt, d, c, mx, primes
from decimal import Decimal as D

l=I(23500); x=log(l)
T=I('.4325');a=I('.1082'); b=I('.324');k=I('3.26');theta=I('.698')
W=exp(k*x);H=exp(theta*x)
hbig=lambda s:exp(-s)*3/s
hsmall=lambda s:exp(-mx(2,s))
Delt=I('9.65')*(k+theta)
assert sum(log(p) for p in primes(1500)).hi<1462
assert ((I('.5')-T)*l-1462-(k+8)*x).lo>0
assert ((I('.5')-T)*l-(k+8)).lo>0
assert (l/(Delt*x)*(1-1/x)).lo>(k+4).hi
assert (k+2-25*theta+1/x).hi<0
assert H.lo>1001
scale=4*l*l*W
low=scale*20*(k*x/2+2)*exp(-l/(Delt*x))
high=scale*8*((k+theta)*x/2+I('2.5'))*I('1.001')*(2/d)**24*exp(-25*theta*x)
triv=scale*2*exp(-l)
# Initial A principal adjustment uses E_small / W, so bound by 1/W times the full component
sml=(1+1/W)*(low+high+triv)
mid=800*sqrt(x)*exp((I('2.5')-k)*x)
big=4000*exp(-I('1.5')*x)
tiny=8*l**3*exp(-l/2)
rem=sml+mid+big+tiny
assert (l-100*x-log(1000)).lo>0 and (3/x+37/log(2)).hi<55 and (51+log(2)/x).hi<52
assert (18*l**3*exp(-l/12)).hi<1 and (5*(1+l/2)*l**7*exp(-l/6)).hi<1
s=T/a
assert s.lo>3 and s.hi<4 and ((T-b)/a).lo>1 and ((T-a)/a).hi<3
bottom=2*exp(I('.577'))*log(s-1)/s-I(134)/300*exp(2)*hbig(s)
def upper_s(s):
    return 2*exp(I('.578'))/s+I(133)/300*exp(2)*hsmall(s)
uppermax=upper_s((T-b)/a)
loss=I(0); weighted=I(0)
for j in range(512):
    A=a+(b-a)*j/512; B=a+(b-a)*(j+1)/512
    # Prime harmonic sum includes a possible prime at either endpoint.
    reciprocal=log(B/A)+1/(2*(A*l)**2)+1/(2*(B*l)**2)+3/(exp(A*l)-2)
    loss+=(b-A)*upper_s((T-B)/a)*reciprocal/(4*b-1)
    weighted+=(b-A)*reciprocal/(4*b-1)
main=(bottom-loss)*I('.26')*I('.60')/a
lam=(b-a)/(4*b-1)
# Coefficient / L for combined error in all weighted divisibility counts, including optional q-character cost 2.
fac=2*mx(I('.90')*a*(1+weighted)+loss/l, lam*(I('.90')*a+uppermax/l))
assert (exp(I('.578'))/2*(1+1/(a*l)**2)).hi<D('.90')
assert fac.hi<D('.3')
deletion=(l*l+2*l)*exp(-l/2)+l*l/(exp(a*l)-1)
assert c.lo>D('.31') and W.lo>400000
assert (a*l-7*x).lo>0 and (a*l-7).lo>0 and k.hi<7
assert (l/2-(7+8)*x).lo>0 and (l/2-15).lo>0
assert (I(8)/9*I('.99')*c-2*(low+high+triv)/scale-l*exp(-l)).lo>D('.26')
assert (I('.0024')*(1-exp(-2))/(c*(l-2-24*d))).hi<D('.01')
assert (I('2.08')/l**2+(2*l+4)*l**49*exp(-l)).hi<1
assert (4*(1+l**14)*l**51*exp(-l)).hi<1
assert (55*x*(l**28+l**14)*exp(-l)).hi<1
assert (exp(1+I('.578')/(51*x))*(I('3.252')*52/14+I('6.505'))).hi<53
assert (1/sqrt(l*x)+sqrt(I('2.08')*56*2*53)).hi<200
assert 98*(4+4+1+1)<1000
assert (I(2)*mx(1,mx(lam*uppermax,lam))).hi<10
assert (I(30)/l**8 + I('.15')*exp(-l/2)).hi<8
# Independently verify the short prime dimension products inherited from linear-tail.md.
eps=1/I(300)
v_core=I(1)
for p in primes(1500):
    v_core*=1-1/I(p-1)
assert (v_core*log(1500)/(1+eps)).lo>D('.60')
sum_logs=I(0)
prefix=None
excess=I(0)
for p in primes(7350):
    if p<=1500:
        continue
    pt=log(log(p))
    prefix=pt-sum_logs if prefix is None else mx(prefix,pt-sum_logs)
    sum_logs+=log(I(p-1)/(p-2))
    excess=mx(excess,prefix+sum_logs-pt)
assert (exp(excess)-1).hi<D('.002902')
assert exp(I('8.9')).hi<7350
assert (exp(I('.001614')+I('.000002964')/log(1500)+2/I(1498))-1).hi<eps.lo
assert rem.hi<D('2') and main.lo>D('.09')
score=main-fac*rem/4-deletion-tiny
assert score.lo>0
print('small:',sml.hi,'middle:',mid.hi,'large:',big.hi,'4 L^2 (E+v) with a negligible extra allowance, upper:',rem.hi)
print('lower dimensionless:',bottom.lo,'detector loss:',loss.hi,'coefficient of R:',fac.hi)
print('main/(N/L):',main.lo,'survivor/(N/L) lower:',score.lo)
