%% 
function y=cfun2(v,a1,a2,a3,a4,k)
C1=[a1(1:(k+1)/2);a1((k+1)/2:-1:2)];
 lam1=fft(C1);
 C2=[a2(1:(k+1)/2);a2((k+1)/2:-1:2)];
 lam2=fft(C2);
 C3=[a3(1:(k+1)/2);a3((k+1)/2:-1:2)];
 lam3=fft(C3);
 C4=[a4(1:(k+1)/2);a4((k+1)/2:-1:2)];
 lam4=fft(C4);
 
n=length(lam1);
v=reshape(v,n,n);
v=fft2(v);
Lam=kron(ones(n,1),lam1)+kron(lam2,ones(n,1))-kron(lam4,lam3);
Lam=1./Lam;
Lamda=reshape(Lam,n,n);
Lamda=Lamda.*v;
Lamda=ifft2(Lamda);
y=Lamda(:);
end