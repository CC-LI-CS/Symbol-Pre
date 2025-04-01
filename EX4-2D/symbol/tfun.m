function y=tfun(v,c1,c2,c3,c4,k)

c1=tau_lambda(c1,k);
c2=tau_lambda(c2,k);
c3=tau_lambda(c3,k);
c4=tau_lambda(c4,k);


n=length(c1);
v=reshape(v,n,n);

V=sqrt(2/(n+1))*dst(v);
V=sqrt(2/(n+1))*dst(V');
V=V';
Lam=kron(ones(n,1),c1)+kron(c2,ones(n,1))-kron(c4,c3);
Lam=1./Lam;
Lam=reshape(Lam,n,n);
Y=Lam.*V;

Y=sqrt(2/(n+1))*dst(Y);
Y=Y';
Y=sqrt(2/(n+1))*dst(Y);
Y=Y';
y=Y(:);
end