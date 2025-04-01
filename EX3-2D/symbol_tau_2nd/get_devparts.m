function [xpdevpart,xmdevpart,ypdevpart,ymdevpart]=get_devparts(p,xgrid,ygrid,xdevord,ydevord,rend)
xpart=(xgrid.^p).*((rend-xgrid).^p);
ypart=(ygrid.^p).*((rend-ygrid).^p);
[Ldev,Rdev]=LR_dev(p,xdevord,xgrid,rend);
xpdevpart=kron(ypart,Ldev);
xmdevpart=kron(ypart,Rdev);

[Ldev,Rdev]=LR_dev(p,ydevord,ygrid,rend);
ypdevpart=kron(Ldev,xpart);
ymdevpart=kron(Rdev,xpart);
end