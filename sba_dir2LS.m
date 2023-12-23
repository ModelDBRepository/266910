function [L,pos]=sba_dir2LS(N,m,gamma,rho,delta,dims,noise_factor)

G=zeros(N,N,'logical');
pos=zeros(N,3);
L=[];
T=(rand(m,m)<rho).*(ones(m,m)-eye(m));
G(1:m,1:m)=T;
for i=1:m
    for j=1:m
        if T(i,j)==1
            L=[L;i,j];
        end
    end
end
targets=[1:m];
dout=sum(T');
din=sum(T);
source=m+1;
D = cumsum(gamma);

pos(1:m,:)=rand(m,3).*dims;

DG=digraph(G);
h=[0];
for i=2:m
    h=[h length(shortestpath(DG,1,i))];
end

wb=waitbar(0, 'inizializzato...');
while source<N
    cc=0;
    lt=length(targets);
    for ii=1:lt
        i=targets(ii);
        G(i,source)=1;
        L=[L;i,source];
        dout(i)=dout(i)+1;
        cc=cc+1;
    end
    dout=[dout,0];
    din=[din, cc];
    DG=digraph(G);
    pos(source,:)=rand(1,3).*dims;
    h=[h length(shortestpath(DG,1,source-1))];
    de=sum((ones(source-1,1)*pos(source,:)-pos(1:(source-1),:)).^2,2);
    %C=delta*(de+200*rand(source-1,1)*noise_factor)+h;
    C=delta*(de+200*rand(source-1,1)*noise_factor)+h(1,1:(source-1))';
    [val,idxs]=sort(C);
    mt=source;
    while mt>source-1
        r=rand();
        mt=find(r<D,1,'first')-1;
    end
    targets=idxs(1:mt);
    source=source+1;
    waitbar(source/N,wb,['Eseguendo ' num2str(round(source/N*100,2)) '%...']);
end
close(wb);