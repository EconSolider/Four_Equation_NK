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
bet=0.99;   %discount factor
kapp=1-1/40; %coupon decay rate
betb=0.99;   %discount factor
ep=11;      %elasticity of substitution of intermediate goods
rhoQE=0.8;
rhor=0;
ph=0.75;     %calvo price
XFI=0.046;   %transfer to bank
rhoA=0.8;
rhothet=0.8;
phipi=1.5;   %Taylor rules
phiy=0.25;
QEs=0.1;     %steady state of QE
rhoZ=0;

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
%19 ZLB
[mcp='Rs>1']
log(Rs)=rhor*log(Rs(-1))+(1-rhor)*(log(steady_state(Rs))+phipi*log(pi/steady_state(pi))
         +phiy*log(Y/steady_state(Y)))-epsr;
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

steady;
check;

shocks;
var epsZ;
periods 1:6;
values 0.2;
 
 var epsQE;
 periods 1;
 values 0.1;
end;

perfect_foresight_setup(periods=200);
perfect_foresight_solver(lmmcp);

for i=1:length(M_.endo_names)
assignin('base', M_.endo_names{i}, oo_.endo_simul(i,:));
end

leg=30;
start=2;

subplot(2,3,1);
plot(Rs(1:leg));
title('Rs');
xlim([start,leg]);

subplot(2,3,2);
plot(QE(1:leg));
title('QE');
xlim([start,leg]);

subplot(2,3,3);
plot(C(1:leg));
title('C');
xlim([start,leg]);

subplot(2,3,4);
plot(Y(1:leg));
title('Y');
xlim([start,leg]);

subplot(2,3,5);
plot(pi(1:leg));
title('\pi');
xlim([start,leg]);

subplot(2,3,6);
plot(Z(1:leg));
title('Z');
xlim([start,leg]);