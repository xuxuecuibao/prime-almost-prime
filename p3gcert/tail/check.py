"""GRH sieve on log N >= 900. Endpoints and scalar envelopes."""
from math import comb, factorial
from decimal import Decimal, ROUND_CEILING
from interval_math import I, log, exp, mx, c, primes

a = I('.1082')
b = I('.324')
T = I('.4325')
s = T/a
odd = primes(1500)
vcore = I(1)
for p in odd:
    vcore *= 1-1/I(p-1)
assert vcore.lo > Decimal('.08')
assert c.lo > Decimal('.31') and c.hi < Decimal('.321')
assert sum(log(p) for p in odd).hi < 1462
assert (vcore*log(1500)/(1+1/I(300))).lo > Decimal('.60')

# Dimension bound for the presieved set of primes, as in the P3 tail proof.
prefix=None
total=I(0)
excess=I(0)
for prime in primes(7350):
    if prime<=1500:
        continue
    position=log(log(prime))
    prefix=position-total if prefix is None else mx(prefix,position-total)
    total+=log(I(prime-1)/(prime-2))
    excess=mx(excess,prefix+total-position)
assert (exp(excess)-1).hi < Decimal('.002902')
assert exp(I('8.9')).hi < 7350
assert (exp(I('.001614')+I('.000002964')/log(1500)+2/I(1498))-1).hi < (1/I(300)).lo

def upper_display(value):
    return value.quantize(Decimal('.000000000001'),rounding=ROUND_CEILING)

def upper(s):
    return 2*exp(I('.578'))/s + I(133)/300*exp(2-mx(2,s))
low = 2*exp(I('.577'))*log(s-1)/s-I(134)/300*exp(2-s)*3/s
assert upper((T-b)/a).hi<5

def run(left, right, B, j):
    L=I(left)
    weighted=I(0)
    loss=I(0)
    for i in range(512):
        u=a+(b-a)*i/512
        v=a+(b-a)*(i+1)/512
        H=log(v/u)+(1/(u*L)**2+1/(v*L)**2)/2+3/(exp(u*L)-2)
        w=(b-u)/(4*b-1)
        weighted+=w*H
        loss+=w*upper((T-v)/a)*H
    assert weighted.hi < Decimal('.49')
    assert (I('.90')*a*(1+weighted)+5*weighted/L).hi < Decimal('.15')
    assert (I('.9')-exp(I('.578'))/2*(1+1/(a*L)**2)).lo>0
    if B is None:
        M=2**len(odd)
        eps=I(0)
        max_core=I(1462)
    else:
        small=[p for p in odd if p<=B]
        tail=[p for p in odd if p>B]
        v=I(1)
        for p in small:
            v *= 1-1/I(p-1)
        eps=v*sum(1/I(p-1) for p in tail)**j/factorial(j)
        M=2**len(small)*sum(comb(len(tail),k) for k in range(j+1))
        max_core=sum(log(p) for p in small)+j*log(1500)
    assert (max_core/L+T).hi<1
    # Each reduced class costs <=12 L sqrt(N); after anchor|A| and at most 3D weights.
    e_char=72*M*L**2*exp((T-I('.5'))*L)
    e_density=I('.1')*eps*I(right if right else 0)**2
    assert (I(5)*exp((b-T)*L)).hi < 1
    shell=eps*I('.33')+12*M*L*exp(-L/2)
    assert (I('.99')*c-shell/vcore).lo>Decimal('.30')
    assert (c*I('1.01')).hi<Decimal('.33')
    assert (I('.0024')*(1-exp(-2))/(c*(L-3))).hi<Decimal('.01')
    assert (42/L).hi<1
    main=(low-loss)*I('.30')*I('.60')/a
    dels=(L*L+2*L)*exp(-L/2) + L**2/(exp(a*L)-1)
    good=main-e_char-e_density-dels
    assert good.lo>Decimal('.04')
    print(left,right, 'main>',str(main.lo)[:13], 'residue error<',upper_display(e_char.hi),
          'density error<',upper_display(e_density.hi), 'good>',str(good.lo)[:14])
run(900,1000,31,8)
run(1000,3000,23,10)
run(3000,None,None,None)
