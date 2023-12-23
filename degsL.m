function [indeg,outdeg]=degsL(L,nn)
% N=size(L,1);
% dv=zeros(1,N);
% indeg=zeros(1,nn);
% outdeg=zeros(1,nn);

% f = waitbar(0,'Please wait...');
% for i=1:N
%     indeg(L(i,2))=indeg(L(i,2))+1;
%     outdeg(L(i,1))=outdeg(L(i,1))+1;
%     waitbar(i/N,f,'Please wait...');
%     disp(i/N)
% end
% f.close()
Lu=unique(L,'rows');
GID=double(unique(Lu(:)));
indeg=histcounts(Lu(:,2),[-0.5,GID'+0.5]);
outdeg=histcounts(Lu(:,1),[-0.5,GID'+0.5]);