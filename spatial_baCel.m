function [Lt,post,ta,blk_id]=spatial_baCel(N,delta,din,dout,Ek,noise_factor,L,pup,pdw)

%close all;
figure(13); hold on;
rho=0.5; 

blk_id=[];
figure(10); hold on; title('gamma');
hdi=histogram(din,'Normalization','pdf');
p_col=[0.98,0.96,0.6;0.98,0.78,0.95;0.96,0.55,0.34;0.77,0.87,0.60;0.4,0.75,0.55;0.56,0.69,0.85];
p_col=p_col(6:-1:1,:);

m=10;
sd=1;
dims=[500 500 2000];
NClayers=[0.5,0.5];

tarr=linspace(0,1,sd+1);
Offx=(ones(sd,sd).*tarr(1:end-1))*dims(1); Offy=(ones(sd,sd).*tarr(1:end-1))'*dims(2);
cs=cumsum(NClayers);
Offz=([0 cs(1:end-1)])*dims(3);
ns=length(Offz);
mm=mean(din) 

h=hdi;

ed=h.BinEdges;
data_vedi=[0]; data_vhi=[0]; c=2;
for i=1:h.NumBins
    if h.Values(i)>0
        data_vedi(c)=(ed(i)+ed(i+1))*0.5;
        data_vhi(c)=h.Values(i);
        c=c+1;
    end
end

nf=floor(data_vedi(end));
hi2=interp1(data_vedi,data_vhi,[0:nf]);
%gamma=hi2;
plot(data_vedi,data_vhi,'k');

T=hi2.*[0:nf];

Es=sum(T); E=cumsum(T);
s=sum(hi2); S=cumsum(hi2);
V=((Es-E)-[1:nf+1].*(s-S))./(s-S);
figure(10); hold on;
plot(V); plot([1 nf+1],[mm-Ek mm-Ek]);
del=find(V<(mm-Ek),1,'first')

APN=N/(sd^2*ns);
AAN=N-APN;
p=(Ek/AAN)

disp(ns*(sd^2))
pop=zeros(1,ns*(sd^2));

for i=1:ns*(sd^2)
    zd=floor((i-1)/(sd^2))+1;
    pop(i)=round(N/sd^2*NClayers(zd));
end
tc=cumsum(pop);
root=[0 tc(1:end-1)];

for ly=1:ns
    PN=round(N/sd^2*NClayers(ly));
    AN=N-PN;
    c=(mm-Ek-m/PN*(m-1)*rho)*PN/(PN-m)
    
    ni=[1:floor(data_vedi(end))];
    
    % sk=cumsum(keri); lambda2=sk(floor(M/3)); lambda1=1-2*lambda2;
    
    hi1=zeros(1,floor(data_vedi(end)));
    % hi1(1:data_vedi(end))=interp1([0 data_vedi],[0 data_vhi],ni);
    hi1(1:data_vedi(end))=interp1(data_vedi,data_vhi,ni);
    
    
    
    alpha_true=zeros(1,ni(end));
    alpha_true([1:ni(end)-del+1])=hi1([del:ni(end)]);
    alpha_true=alpha_true/sum(alpha_true);
    sa=sum(alpha_true.*[0:ni(end)-1]);
    la=length(alpha_true);
    
    iota=binopdf([0:la-1],m-1,rho);
    gamma{ly}=PN/(PN-m)*alpha_true-m/(PN-m)*iota;
    gamma{ly}=gamma{ly}/sum(gamma{ly});
end



% [G,pos]=sba_dir2(N,m,gamma,rho,delta,[1 1 1]);
% figure(2); hold on; title('distances');
% dvs=dist_vec(G,pos);

Lt=[];

post=[];
ta=zeros(1,ns*(sd^2));
%pool_obj=openPool('local',4);

Ltt={};
parfor i=1:ns*(sd^2)
    Lttw=[];
    for j=1:ns*(sd^2)
        if i==j
            zd=floor((i-1)/(sd^2))+1;
            offp=pop(i);
            xyd=mod(i-1,(sd^2))+1;
            [xd,yd]=ind2sub(sd,xyd);
            tic;
            [Lr,posr]=sba_dir2L(offp,m,gamma{zd},rho,delta,[1/sd 1/sd NClayers(zd)].*dims,noise_factor);
            toc
            tp=[Offx(xd,yd)+posr(:,1),Offy(xd,yd)+posr(:,2),Offz(zd)+posr(:,3)];
            post=[post; tp];
            blk_id=[blk_id,ones(1,size(tp,1))*i];
            Lr(:,1)=Lr(:,1)+root(i);
            Lr(:,2)=Lr(:,2)+root(j);
            figure(5); hold on; title('positions');
            plot3(tp(:,1),tp(:,2),tp(:,3),'x','color',p_col(zd,:));
        else
            %Gr=(rand(N/(ns*4),N/(ns*4))<0.0005);
            %Lr=[];
%             for ri=1:pop(i)
%                 for rj=1:pop(j)
%                     if rand()<p
%                         Lr=[Lr;ri+root(i),rj+root(j)];
%                     end
%                 end
%             end
            rpi=floor(pop(i)/L); rpj=floor(pop(j)/L)
            MT=(rand(rpi,rpj)<p);
            MT2=zeros(pop(i),pop(j));
            MTpos=find(MT==1); MTzero=find(MT==0);
            for ti=1:L
                for tj=1:L
                    TMT=MT;
                    TMT(MTpos)=(rand(1,length(MTpos))<pup);
                    TMT(MTzero)=(rand(1,length(MTzero))<pdw);
                    MT2(L*(0:(rpi-1))+ti,L*(0:(rpj-1))+tj)=TMT;
                end
            end
            permi=randperm(pop(i),pop(i)); permj=randperm(pop(j),pop(j));
            [Ltr,Ltc]=find(MT2(permi,permj));
            Lr=[Ltr+root(i),Ltc+root(j)];
            %sum(sum(Gr))
        end
        %Gtt=[Gtt,Gr];
        Lttw=[Lttw;Lr];
        %Lra{j}=Lr;
        %disp(j)
    end
    disp(i)
    Ltt{i}=Lttw;
%     disp('recuperando parti!')
%     for j=1:ns*(sd^2)
%         Ltt=[Ltt;Lra{j}];
%     end
    %Gt=[Gt;Gtt];
    disp(['iterazione: ' num2str(i)]);
end


for i=1:ns*(sd^2)
    Lt=[Lt;Ltt{i}];
end


size(Lt)

legend({'sim','data'},'Location','northeast');


[indeg,outdeg]=degsL(Lt,N); 
totdeg=indeg+outdeg;
mi=max(indeg); mo=max(outdeg); mt=max(totdeg);
d=sum(indeg)-sum(outdeg)

% figure(1); hold on;
% subplot(1,2,1); hold on;
% sci=(0.001+indeg/mi)*30;
% title('outdegree scatter plot');
% scatter3(post(:,1),post(:,2),post(:,3),sci,indeg);
% colorbar
% subplot(1,2,2); hold on;
% sco=(0.001+outdeg/mo)*30;
% title('outdegree scatter plot');
% scatter3(post(:,1),post(:,2),post(:,3),sco,outdeg);
% colorbar
% 
% figure(2); hold on;
% sct=(0.001+(totdeg)/mt)*30;
% title('totdegree scatter plot');
% scatter3(post(:,1),post(:,2),post(:,3),sct,totdeg);
% colorbar
% axis([ 0 dims(1) 0 dims(2) 0 dims(3)])

nsud=50;
figure(3); hold on;
subplot(1,2,1); hold on; title('indegree');
hi=histogram(indeg,nsud,'Normalization','pdf');
hid=histogram(din,nsud,'Normalization','pdf');
legend({'sim','data'},'Location','northeast');

subplot(1,2,2); hold on; title('outdegree');
ho=histogram(outdeg,nsud,'Normalization','pdf');
hout=histogram(dout,nsud,'Normalization','pdf');
legend({'sim','data'},'Location','northeast');

