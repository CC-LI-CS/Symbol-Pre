function [Ldev,Rdev]=LR_dev(p,ord,grid,rend)
rms=rend-grid;
sppart=grid.^(p-ord);
rmsppart=rms.^(p-ord);
weights=nchoosek(p,0)*(rend^p)*gamma(p+1)/gamma(p+1-ord);
Ldev=weights*sppart;
Rdev=weights*rmsppart;
for kk=1:p
    weights=nchoosek(p,kk)*(rend^(p-kk))*((-1)^kk)*gamma(kk+p+1)/gamma(kk+p+1-ord);
    sppart=grid.*sppart;
    rmsppart=rms.*rmsppart;
    Ldev=Ldev+weights*sppart;
    Rdev=Rdev+weights*rmsppart;
end
end