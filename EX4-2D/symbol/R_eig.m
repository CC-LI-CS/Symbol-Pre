% Riesz preconditioner
function [RES1,RES2]=R_eig(afa_1,afa_2,k)

g1=wfun(afa_1,k);
g2=wfun(afa_2,k);
% g3=-wfun(afa_3,k);
% g4=-wfun(afa_4,k);

RES1=tau_lambda(g1,k);
% RES1=RES1;
RES2=tau_lambda(g2,k);
% RES2=RES2;
% RES3=tau_lambda(g3,k);
% RES4=tau_lambda(g4,k);