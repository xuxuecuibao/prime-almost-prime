"""Four-factor certificate for 91 <= log N <= 115."""
from parameter_tables import SUPPORT_RECIPES, iter_parameter_rows
from decimal import Decimal,ROUND_FLOOR
from p4_fd_verify import I,exp,log,mn,mx,branch,local,P,K,H,d,c
D=Decimal
def floor_decimal(x):
 return x.to_integral_value(rounding=ROUND_FLOOR)
def check(row,recipes):
 l0,u0,B,b,j,mh,md,mf=row;l,u=I(l0),I(u0);h=I(mh)/1000;dep=I(md)/1000;fr=I(mf)/100
 assert 91<=l0<u0<=115 and 2<b<=B<=20000 and j in (2,4)
 assert (5*h-1).lo>D('.5') and h.hi<D('.5') and fr.lo>=0
 y,v,delta,T,square,q,rho=local(B,b,j)
 eta=K/c*(20*(log(K)/2+2)*exp(-l/2)+(log(K*(H+1))/2+2)*8/(H-1)*(2/(d*H))**24+2*exp(-l))
 xs=[log(K)];es=[eta];ac=I(0);old=1/I(K)
 for i,ch in enumerate(recipes,1):
  a,t,s,nu=map(D,ch)
  assert a>0 and t>1 and s>0 and nu>1
  qq=log(K)+i*I('.12');z=exp(-qq); A=branch(l,u,qq,ch);ac+=A*(old-z);old=z
  xs.append(qq);es.append(eta+(ac+A*z)/c)
 for i in range(len(es)-2,-1,-1):es[i]=mn(es[i],es[i+1])
 def err(qx):
  from bisect import bisect_left
  k=bisect_left([z.lo for z in xs],qx.hi)
  assert k<len(es);return es[k]
 m1=(1-exp(-2))/c*I('.0024')/(l-2-24*d)
 s0=(1-m1)*(1-T[j+2]/v)-err(q+(j+1)*y)/delta-I('1e-7')
 assert s0.lo>0
 den=c*delta*s0
 # For p²<=K restrict to the largest core primorial P0 giving p²P0<=K.
 # Apply Platt to the reduced CRT classes, with weighted relative error <=m1+eta+1e-7.
 # Assign small sieve primes dividing N an arbitrary forbidden reduced class.
 sq0=I(0);sleft=square
 for p in P:
  if p<=B:continue
  if p*p>K:break
  Q=1;V=I(1)
  for s in P:
   if s>B or Q*p*p*s>K:break
   Q*=s;V*=1-1/I(s-1)
  g=1/I(p*(p-1));sq0+=(1+m1+eta+I('1e-7'))*g*V/(s0*delta);sleft-=g
 assert sleft.lo>0
 squares=sq0+2*sleft/(I('.6')*den)+2*u/((exp(l/5)-1)*den)
 assert (K*20*u*exp(-l)/c).hi<D('1e-7')
 last=I(0);maxQ=I(0);S=I(0);W=I(0)
 for n in range(j,-1,-2):
  support=dep-fr*n*y
  ts=(support-q-n*y)/u
  t=I(min(D(mh)/1000,floor_decimal(ts.lo*10**8)/10**8))
  if t.hi<=max(last.lo,(y/l).lo):continue
  assert t.lo>max(last.hi,(y/l).hi)
  if last.hi==0:F=h*(log(t*u/y)+rho)-t+y/l
  else:
   w=last*l
   assert w.lo>log(286).hi, (l0,n,w.lo)
   F=h*log(t/last)-t+last+(h-last)*(1/w**2+1/(exp(w)-1))
  assert F.lo>0
  sT=local(B,b,n)[3];gap=sT[n+1]/v
  ss=q+(n+1)*y;pp=(I(B)**(n+1))*(u*exp(-l/2)+(b+1)*u*exp(-l))/den
  assert ss.hi<l.lo
  R=mn(((1+m1)*(1+gap)+err(q+n*y)/delta+I('1e-7'))/s0,1+2*gap/((1-ss/l)*c*s0)+pp+I('1e-7'))
  S+=F*R
  a=mx(y/u,last)
  for sub in range(1,9):
   sx=a+(t-a)*sub/8
   W+=(t-a)/8*err(mx(maxQ,q+n*y+sx*u))
  maxQ=mx(maxQ,q+n*y+t*u);last=t
  if t.lo==t.hi==h.lo or t.lo>=D(mh)/1000:break
 assert abs(last.lo-D(mh)/1000)<=D('1.000001e-8')
 if last.lo < D(mh)/1000:raise ArithmeticError("clip h exactly")
 indu=(b+j+2)*u*(exp(q+(j+1)*y-l)+2*exp(maxQ-l))/den
 del0=u*(exp(-l/2)+(B+3+u)*exp(-l))/den
 assert indu.hi<D('1e-9') and del0.hi<D('1e-7')
 bad=(S+W/(s0*delta))/(5*h-1)+squares+I('.000001')
 print(l0,u0,'candidate',format(s0.lo,'.8f'),'excluded',format(bad.hi,'.8f'),flush=True)
 assert bad.hi<D('.992878')
 return bad.hi
def main():
    for row, choices in iter_parameter_rows(SUPPORT_RECIPES):
        check(row, choices)

if __name__=='__main__':main()
