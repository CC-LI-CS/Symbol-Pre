function r1=fouriercoefficient(k,alpha,n)
v=(1:k-1)';
b=sin(pi/2*v)./v;
b=1/pi*b;
h=pi/n/2;
 
c1=v.^(1+alpha);
c1=(1/pi)./c1;

a1=firstcolumn(n,h,alpha,k);

a0=pi^alpha/2^(alpha+1)/(alpha+1)+0.5;

a1=c1.*a1-b;
r1=[a0;a1]; % first column of toeplitz in x direction
end