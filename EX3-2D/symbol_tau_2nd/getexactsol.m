function vout=getexactsol(xgrid,ygrid,p,rend)
vout=kron((ygrid.^p).*(rend-ygrid).^p,(xgrid.^p).*(rend-xgrid).^p);

end