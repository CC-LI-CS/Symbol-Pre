% compute  Toelitz matrix system kron(I,A)+kron(B,I)+kron(C,D)
clc;
clear;
n=1001;   % 
for m=7
k=2^m-1;
afa_1=1.01;
afa_2=1.01;
afa_3=1;
afa_4=1;    %1/k^(afa_2)*
a1=fouriercoefficient(k,afa_1,n);
a2=fouriercoefficient(k,afa_2,n);
a3=fouriercoefficient(k,afa_3,n);
a4=fouriercoefficient(k,afa_4,n);

Time_start=cputime;

%% 2D toeplitz matrix multiply vector
c1=[a1;0;a1(k:-1:2)];
c1=fft(c1); % the eigenvalues of Toeplitz matrix in X direction

c2=[a2;0;a2(end:-1:2)];
c2=fft(c2); 

c3=[a3;0;a3(k:-1:2)];
c3=fft(c3); % 

c4=[a4;0;a4(k:-1:2)];
c4=fft(c4); % 


%%
 s=rng;
u_e=rand(k^2,1);
rng(s);
%u0=zeros(k^2,1);
u0 = ones(k^2,1)/sqrt(k^2);
size(u_e)
%%
   %[u,flag,reles,iter]=pcg(@(z)afun(c1,c2,c3,c4,z),afun(c1,c2,c3,c4,u_e),1e-8,1000,[],[],u0);% N_pre
  tic
   [u,flag,reles,iter]=pcg(@(z)afun(c1,c2,c3,c4,z),afun(c1,c2,c3,c4,u_e),1e-8,1000,@(p)Rfun(p,afa_1,afa_2,k),[],u0); % R_pre
    iter3=iter
    toc
    tic
    [u,flag,reles,iter]=pcg(@(z)afun(c1,c2,c3,c4,z),afun(c1,c2,c3,c4,u_e),1e-8,1000,@(p)tfun(p,a1,a2,a3,a4,k),[],u0);%tau_pre
    iter2=iter
toc
tic
 [u,flag,reles,iter]=pcg(@(z)afun(c1,c2,c3,c4,z),afun(c1,c2,c3,c4,u_e),1e-8,1000,@(p)cfun2(p,a1,a2,a3,a4,k),[],u0);%C_pre
iter1=iter
toc
%err(m-9)=norm(u-u_e,inf)
%it(m-9)=iter
%CPU_time(m-9)=cputime-Time_start
%it

end

