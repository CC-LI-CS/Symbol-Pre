clc;
clear;
d1p=2;d1m=35;
d2p=1;d2m=20;
rend=2;
xl=0;xr=rend;
yl=0;yr=rend;
p=2;

levels=[7 9];
alpha1s=[1.1];
alpha2s=[1.1];
T=1;


for iii=1:length(alpha1s)
    alpha1=alpha1s(iii);
    for hhh=1:length(alpha2s)
    alpha2=alpha2s(hhh);
    disp('(alpha1,alpha2)=')
    disp([alpha1,alpha2])
    iternum=[];
    t=[];
    err=[];
    DoF=[];
    for jjj=1:length(levels)
        level=levels(jjj);
        M=2^level;
        N=2^level-1;
N1=N;
N2=N;
tau=T/M;
DoF(jjj)=N1*N2;
h1=(xr-xl)/(N1+1);
h2=(yr-yl)/(N2+1);
eta1p=d1p;eta1m=d1m;eta2p=d2p*(h1^alpha1)/(h2^alpha2);eta2m=d2m*(h1^alpha1)/(h2^alpha2);
nu=2*(h1^alpha1)/tau;

xgrid=(h1:h1:xr-h1).';
ygrid=(h2:h2:yr-h2).';
tgrid=0:tau:T;
% tpgrid=tgrid+0.5*tau;

tic;

Txseries=get_Tseries(alpha1,N1);%1ST order SGD
Tyseries=get_Tseries(alpha2,N2);

xMatfcol=eta1p*Txseries(2:end);
xMatfcol(1:2)=xMatfcol(1:2)+eta1m*[Txseries(2);Txseries(1)];%column of nonsysmmetric toeplitz
xMatfrrow=eta1m*Txseries(3:end); xMatfrrow(1)=xMatfrrow(1)+eta1p*Txseries(1);%row of nonsysmmetric toeplitz
augTxeigs=get_augmented_toep_eig(xMatfcol,xMatfrrow);%the eigenvalue of nonsysmmetric toeplitz
stxfcol=(eta1p+eta1m)*[Txseries(2);0.5*(Txseries(1)+Txseries(3));0.5*Txseries(4:end)];%(G+G^T)/2
tauxeigs=get_taumateig(stxfcol);

yMatfcol=eta2p*Tyseries(2:end);
yMatfcol(1:2)=yMatfcol(1:2)+eta2m*[Tyseries(2);Tyseries(1)];
yMatfrrow=eta2m*Tyseries(3:end); yMatfrrow(1)=yMatfrrow(1)+eta2p*Tyseries(1);
augTyeigs=get_augmented_toep_eig(yMatfcol,yMatfrrow);
styfcol=(eta2p+eta2m)*[Tyseries(2);0.5*(Tyseries(1)+Tyseries(3));0.5*Tyseries(4:end)];
tauyeigs=get_taumateig(styfcol);

s1fac=sqrt(2/(N1+1));s2fac=sqrt(2/(N2+1));
scasqrtaueigs=nu+(kron(ones(N2,1),tauxeigs)+kron(tauyeigs,ones(N1,1)));
scasqrtaueigs=reshape((s1fac*s2fac)./sqrt(scasqrtaueigs),N2,N1);%%%%%%
[xpdevpart,xmdevpart,ypdevpart,ymdevpart]=get_devparts(p,xgrid,ygrid,alpha1,alpha2,rend);
uexact=getexactsol(xgrid,ygrid,p,rend);

f1=sourceterm(d1p,d1m,d2p,d2m,xpdevpart,xmdevpart,ypdevpart,ymdevpart);
f2=2*(h1^alpha1)*exp(tau/2)*(uexact+f1);
f3=rhsminresmv(nu,augTxeigs,augTyeigs,scasqrtaueigs,N1,N2,uexact);
f=f2+f3;
%% 
f=f(:);
flip(f);
r0=norm(f,2);
restart=20;
maxit=ceil((N1*N2)/restart);
tol=1e-8;

%% 
maxit=ceil((N1*N2)/restart);
 [uapp,flag,~,iter,resvec] =minres(@(x)minresmv(nu,augTxeigs,augTyeigs,scasqrtaueigs,N1,N2,x),flip(f),tol,1000,@(pv)minresinpv(scasqrtaueigs,N1,N2,pv)); % minres 2nd
iternum(jjj)=iter;
t(jjj)=toc;

err(jjj)=norm(uapp(:)-exp(tau)*uexact(:),inf);
    end
    levels
DoF
iternum
t
err
end
end
