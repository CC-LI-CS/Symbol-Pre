%% Toeplitz matrix kronecker identity matrix multiply vector
function y=Ty_K_I_V(lamda,u) %lamda is the eigenvalue of [a;0;a(n1:-1:2)]
n=length(u);
n1=sqrt(n);
U=reshape(u,n1,n1);
U=U';
U=[U;zeros(n1,n1)];
U=fft(U);
Lam=repmat(lamda,1,n1);
Lam=Lam.*U;
Lam=ifft(Lam);
Y=Lam(1:n1,:);
Y=Y';
y=Y(:);
end