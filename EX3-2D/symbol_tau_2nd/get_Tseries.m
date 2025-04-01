%% second order
function toe=get_Tseries(alpha,msize)
% alpha=1.2
% msize=5

%% second order
toe=zeros(msize+1,1);toe(1)=1;
w=zeros(msize+1,1);w(1)=alpha/2;
for k=2:msize+1
   toe(k)=(1-(alpha+1)./(k-1))*toe(k-1); 
    w(k)=(alpha/2)*toe(k)+(2-alpha)/2*toe(k-1);
end
toe=-w;
end


% %% first order
% function toe=get_Tseries(alpha,msize)
% toe=zeros(msize+1,1);toe(1)=-1;
% for k=2:msize+1
%    toe(k)=(1-(alpha+1)./(k-1))*toe(k-1); 
% end
% end