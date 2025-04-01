function vout=gmres_afun(D1p,D2p,D1m,D2m,augTxpeigs,augTypeigs,augTxmeigs,augTymeigs,N1,N2,vin)
% vin=dstn(dstn(scasqrtaueigs.*reshape(vin,N2,N1),1),2);
vin=reshape(vin,N2,N1);
extxlen=length(augTxpeigs(:));
extylen=length(augTypeigs(:));
temp=[vin;zeros(extxlen-N1,N2)];
partx=ifft(fft(temp,[],1).*repmat(augTxpeigs(:),1,N2),[],1);
partx=real(partx(1:N1,:));
t=(partx);
tp1=D1p.*t(:);
%%%%%


partx=ifft(fft(temp,[],1).*repmat(augTxmeigs(:),1,N2),[],1);
partx=real(partx(1:N1,:));
t=(partx);
tm1=D1m.*t(:);

%% 
temp=[vin,zeros(N1,extylen-N2)];
vout=ifft(fft(temp,[],2).*repmat(augTypeigs(:).',N1,1),[],2);
vout=real(vout(:,1:N2));
vout=vout;
tp2=D2p.*vout(:);

vout1=ifft(fft(temp,[],2).*repmat(augTymeigs(:).',N1,1),[],2);
vout1=real(vout1(:,1:N2));
vout1=vout1;
tm2=D2m.*vout1(:);
%% 
% vout=reshape(vout,N1,N2);
% vout=scasqrtaueigs.*dstn(dstn(vout+partx,1),2);
vout=tm1+tm2+tp1+tp2;

end