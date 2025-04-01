%%  identity matrix kronecker Toeplitz matrix multiply vector
function y=I_K_Tx_V(lamda,u) % lamda is the eigenvalue of toeplitz matrix
n=length(u);
n1=sqrt(n);
U=reshape(u,n1,n1);
U=[U;zeros(n1,n1)];
U=fft(U);
% Lam=repmat(lamda,1,n1);
Lam=bsxfun(@times,lamda,U);
% Lam=Lam.*U;
Lam=ifft(Lam);
Y=Lam(1:n1,:);
y=Y(:);
end