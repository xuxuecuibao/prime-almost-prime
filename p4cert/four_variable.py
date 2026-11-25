"""P4 on 212 <= log N <= 460. Requires adjacent fixed_platt_rz.py."""
from decimal import Decimal as D, ROUND_FLOOR
from fixed_platt_rz import K,c,box,H,middle,sieve_numbers,rbox,one

def four_row(l,u,B,b,jb,mh):
    l,u=map(D,(l,u));h=D(mh)/1000;y=D(B).ln();q,v,delta,X,square,rho=sieve_numbers(B,b,jb)
    assert 5*h>1 and X[1]<2
    m1=(1-(-D(2)).exp())/c*D('.0024')/(l-2-rbox*box)
    eta=K*(20*(D(K).ln()/2+2)*(-l/2).exp()+((D(K)*(H+1)).ln()/2+2)*8/(H-1)*(2/(box*H))**rbox+2*(-l).exp())/c
    outer=98*(4/u**D('2.5')+4/l**D('2.5')+18*l.sqrt()*(-l/12).exp()+5*(1+u/2-6*u.ln()-6*l.ln())*l**D('3.5')*(-l/6).exp())/c
    E=eta+middle(l,u,K)/c+outer
    P=(eta if B<K else 0)+middle(l,u,max(B,K))/c+outer
    s0=(1-m1)*(1-X[jb+2]/v)-E/delta-D('1e-7')
    level=l/2-6*l.ln(); assert s0>0 and q+(jb+1)*y<level
    denominator=c*s0*delta
    assert (D(B)**(jb+1)*u*(-l/2).exp())/denominator<D('1e-9')
    assert (4*(b+jb+2)*(u+2)*(-l/2).exp()+2*u/((l/5).exp()-1))/denominator<D('1e-8')
    S=D(0); last=D(0)
    for j in range(jb,-1,-2):
        # Rounding down keeps the progression levels within the proved domain.
        t=min(h,(((level-q-j*y)/l)*10**8).to_integral_value(rounding=ROUND_FLOOR)/10**8)
        if t<=max(last,y/l):continue
        assert q+j*y+t*l<=level and D('.5')-t>6/l
        if last==0:
            F=h*((t*u/y).ln()+rho)-t+y/l
        else:
            a=last*l; F=h*(t/last).ln()-(t-last)+(h-last)*(1/a**2+1/(a.exp()-1))
        assert F>=0
        g=X[j+1]/v
        R=min(((1+m1)*(1+g)+E/delta+D('1e-7'))/s0,1+2*g/((1-(q+(j+1)*y)/l)*c*s0)+D('1e-7'))
        S+=R*F;last=t
        if last==h:break
    assert last==h
    bad=(S+h*P/(s0*delta))/(5*h-1)+2*square/(D('.6')*denominator)+D('1e-6')
    print(int(l),int(u),B,b,jb,mh,f'{s0:.7f}',f'{bad:.7f}',flush=True)
    assert bad<D('.99')

ROWS=[(212,217,700,7,6,320),(217,230,1000,7,8,320),(230,242,1500,7,8,330),(242,266,2500,7,8,340),
      (266,300,2500,7,8,350),(300,340,5000,7,8,370),(340,368,30000,11,8,370),(368,460,30000,11,10,380)]
if __name__=='__main__':
    for params in ROWS:four_row(*params)
