function y=Rfun(v,afa_1,afa_2,k)
%[lambda1,lambda2,lambda3,lambda4]=R_eig(afa_1,afa_2,afa_3,afa_4,k);
[lambda1,lambda2]=R_eig(afa_1,afa_2,k);
n=length(lambda1);
v=reshape(v,n,n);

V=sqrt(2/(n+1))*dst(v);
V1=sqrt(2/(n+1))*dst(V');
V1=V1';
Lam=kron(ones(n,1),lambda1)+kron(lambda2,ones(n,1));
Lam=1./Lam;
Lam=reshape(Lam,n,n);
Y=Lam.*V1;

Y=sqrt(2/(n+1))*dst(Y);
Y1=Y';
Y1=sqrt(2/(n+1))*dst(Y1);
Y1=Y1';
y=Y1(:);
end