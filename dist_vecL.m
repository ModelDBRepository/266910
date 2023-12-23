function dv=dist_vecL(L,posI)


% f = waitbar(0,'Please wait...');
% for i=1:N
%     temp1=posI(L(i,1),:);
%     temp2=posI(L(i,2),:);
%     dv(i)=sqrt(sum((temp1-temp2).^2));
%     waitbar(i/N,f,'Please wait...');
%     disp(100*i/N)
% end
% close(f);

temp1=posI(L(:,1),:);
temp2=posI(L(:,2),:);
d2=(temp1-temp2).^2;
dv=sqrt(sum(d2,2));
%size(dv)


%figure(1); hold on;
hd=histogram(dv,'Normalization','pdf');