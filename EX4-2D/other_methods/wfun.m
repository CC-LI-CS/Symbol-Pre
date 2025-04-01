function w=wfun(alpha,n)
g0=1;
g=zeros(n,1);
g(1)=-alpha;
for j=1:n-1
g(j+1)=(1-(alpha+1)/(j+1))*g(j);
end


% r1=(alpha^2+3*alpha+2)/12;
% r0=(4-alpha^2)/6;
% r_1=(alpha^2-3*alpha+2)/12;
% w0=r1*g0;
% w=zeros(n,1);
% w(1)=r1*g(1)+r0*g0;
% w(2)=r1*g(2)+r0*g(1)+r_1*g0;
% w(3:n)=r1*g(3:n)+r0*g(2:n-1)+r_1*g(1:n-2);
% w(1)=2*w(1);
% w(2)=w0+w(2);
 w=g;
 w(1)=2*g(1);
 w(2)=g0+g(2);

end