var
C lambd L w Rs pi cb b lambdb Cb
Rb Q bFI re s QE thet v omeg Rre
A mc Fa Fb Y pistar bcb Xb Z
;
% Note that Cb is non-labor households' consumption
%           and cb is newly issued bond


varexo epsA epsthet epsQE epsr epsZ;

parameters
sigm ps ch bet kapp betb ep rhoQE rhor ph XFI rhoA rhothet 
phipi phiy QEs rhoZ;
%注意eps是matlab中的固有函数，
%因此用eps代表epsilon会导致稳态问题。

sigm=1;      %inverse elasticity of substitution
ps=0.005;    %labor disutility,！！re-calibarte in steady file！！
ch=1;        %inverse Firsch elasticity
bet=0.995;   %discount factor
kapp=1-1/40; %coupon decay rate
betb=0.99;   %discount factor
ep=11;      %elasticity of substitution of intermediate goods
rhoQE=0.8;
rhor=0.8;
ph=0.75;     %calvo price
XFI=0.046;   %transfer to bank
rhoA=0.8;
rhothet=0.8;
phipi=2.5;   %Taylor rules
phiy=0.25;
QEs=0.1;     %steady state of QE
rhoZ=0.8;

model;
%1
C^(-sigm)=lambd;
%2
ps*L^ch=lambd*w;
%3
lambd=Z*bet*lambd(+1)*Rs/pi(+1);
%4
cb=b-kapp*b(-1)/pi;
%5
lambdb=Cb^(-sigm);
%6
Rb=(1+kapp*Q)/Q(-1);
%7
lambdb=betb*lambdb(+1)*Rb(+1)/pi(+1);
%8
Q*(bFI-kapp*bFI/pi)+re=s+XFI;
%9
Q*bFI=thet*XFI;
%10
lambd(+1)/lambd/pi*(Rb(+1)-Rs)=omeg;
%11
Rre=Rs;
%12
w=A*mc;
%13
Fa=lambd*mc*Y+bet*ph*(Fa(+1))* (pi(+1)^ep);
%14
Fb=lambd*Y+bet*ph*(Fb(+1))* (pi(+1)^(ep-1));
%15
pistar=ep/(ep-1)*pi*Fa/Fb;
%16
re=Q*bcb;
%17
QE=Q*bcb;
%18
log(QE)=rhoQE*log(QE(-1))+(1-rhoQE)*log(steady_state(QE))+epsQE;
%19
[name = 'Nominal Interest Rate', relax='elb']
log(Rs)=(log(steady_state(Rs))+phipi*log(pi/steady_state(pi))
         +phiy*log(Y/steady_state(Y)))-epsr;
[name = 'Nominal Interest Rate', bind='elb']
Rs=steady_state(Rs)-1e-8;
%20
pi^(1-ep)=ph+(1-ph)*pistar^(1-ep);
%21
v*Y=A*L;
%22
v=ph*v(-1)*pi^ep+(1-ph)*pistar^(-ep);
%23
Y=C+Cb;
%24
Xb=(1+kapp*Q)/pi*b(-1);
%25
Cb=Q*b;
%26
b=bFI+bcb;
%27
log(thet)=rhothet*log(thet(-1))
         +(1-rhothet)*log(steady_state(thet))+epsthet;
%28
log(A)=rhoA*log(A(-1))
       +(1-rhoA)*log(steady_state(A))+epsA;
%29 discount factor shock
log(Z)=rhoZ*log(Z(-1))+epsZ;
end;

occbin_constraints;
name 'elb'; 
bind Rs <= steady_state(Rs)-1e-8; %给程序一些容错空间，更容易得到解
relax Rs > steady_state(Rs)+1e-8;
end;

steady;
check;

shocks(surprise,overwrite);
var epsZ;
periods 1;
values 0.05;

% var epsr;
% periods 5;
% values 0.5;

% var epsQE;
% periods 5;
% values 1;
end;

occbin_setup(simul_periods = 20, simul_check_ahead_periods = 400,simul_maxit = 100);
occbin_solver;
occbin_graph Z Rs C Y pi;