% compute Toeplitz block multipply vector
function y=afun(a1,a2,a3,a4,v)
n=length(a1)/2;
y1=I_K_Tx_V(a1,v);
y2=Ty_K_I_V(a2,v);
y3=Ty_K_I_V(a4,v);

U=reshape(y3,n,n);
%U=U';
u=U(:);
y4=I_K_Tx_V(a3,u);
% V=reshape(y4,n,n);
% V=V';
% y4=V(:);
y=y1+y2-y4;
end
