function vout=rhsminresmv(nu,augTxeigs,augTyeigs,scasqrtaueigs,N1,N2,vin)
% vin=dstn(dstn(scasqrtaueigs.*reshape(vin,N2,N1),1),2);
vin=reshape(vin,N2,N1);
extxlen=length(augTxeigs(:));
extylen=length(augTyeigs(:));
temp=[vin;zeros(extxlen-N1,N2)];
partx=ifft(fft(temp,[],1).*repmat(augTxeigs(:),1,N2),[],1);
partx=real(partx(1:N1,:));
temp=[vin,zeros(N1,extylen-N2)];
vout=ifft(fft(temp,[],2).*repmat(augTyeigs(:).',N1,1),[],2);
vout=real(vout(:,1:N2));
vout=nu*vin-(vout+partx);
% pp=size(vout)
% vout=flip(vout);  

%% 
% vout=reshape(vout,N1,N2);
% vout=scasqrtaueigs.*dstn(dstn(vout+partx,1),2);
vout=vout(:);

end