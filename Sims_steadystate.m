function [ys,params,check] = Sims_steadystate(ys,exo,M_,options_)
% Inputs: 
%   - ys        [vector] vector of initial values for the steady state of the endogenous variables
%   - exo       [vector] vector of values for the exogenous variables
%   - M_        [structure] Dynare model structure
%   - options   [structure] Dynare options structure
%
% Output: 
%   - ys        [vector] vector of steady state values for the the endogenous variables
%   - params    [vector] vector of parameter values
%   - check     [scalar] 0 if steady state computation worked and to
%                        1 of not (allows to impose restrictions on parameters)
check=0;
%% Step 1: read out parameters to access them with their name
for ii = 1:M_.param_nbr
  eval([ M_.param_names{ii} ' = M_.params(' int2str(ii) ');']);
end

%% Step 2: Enter model equations here
thet=5;
XFI=0.046;
pi=1;
pistar=1;
A=1;
v=1;
mc=(ep-1)/ep;
L=1;
Rs=1/bet;
Rb=1/betb;
Q=1/(Rb-kapp);
omeg=Rb-Rs;
Rre=Rs;
QE=QEs;
bcb=QE/Q;
re=QE;
bFI=thet*XFI/Q;
b=bcb+bFI;
Cb=Q*b;
Y=1;
C=Y-Cb;
lambd=C^(-sigm);
lambdb=Cb^(-sigm);
s=Q*(1-kapp)*bFI+re-XFI;
Fa=lambd*mc*Y/(1-bet*ph);
Fb=lambd*Y/(1-bet*ph);
w=mc;
Xb=(1+kapp*Q)*b;
ps=lambd*w;
cb=b-kapp*b;

%% Step 3: Update parameters and variables
params=NaN(M_.param_nbr,1);
for iter = 1:M_.param_nbr %update parameters set in the file
  eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ])
end

for ii = 1:M_.orig_endo_nbr %auxiliary variables are set automatically
  eval(['ys(' int2str(ii) ') = ' M_.endo_names{ii} ';']);
end
end