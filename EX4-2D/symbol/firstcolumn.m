function y=firstcolumn(n,h,alpha,k)

x1=(0:n-1)'*h;

x2=(0.5:n-0.5)'*h;
x3=(1:n)'*h;

% y=zeros(k-1,1);
%for j=1:k-1


    j0=pi/2.*(0:k-2);
    x11=bsxfun(@plus,j0,x1);
    x21=bsxfun(@plus,j0,x2);
    x31=bsxfun(@plus,j0,x3);
    f1=cos(x11).*(x11.^alpha);
    f2=cos(x21).*(x21.^alpha);
    f3=cos(x31).*(x31.^alpha);
    f=h/6*(f1+4*f2+f3);
    y=sum(f,1);
%end
y=cumsum(y);
y=y';
end