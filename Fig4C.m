close all
clear all

Iavg=5;


mode='Slice6';
S=300;
Nv=0.1;
[tszSl,tnrSl,dataSl,pathSl]=routine(S,mode,Nv,Iavg);

mode='Prob_Dist';
S=300;
Nv=0.1;
[tszPD,tnrPD,dataPD,pathsPD]=routine(S,mode,Nv,Iavg);

mode='ER';
S=300;
Nv=0.1;
[tszER,tnrER,dataER,pathER]=routine(S,mode,Nv,Iavg);

dbase={dataER,dataPD,dataSl};
slabel={'ER','ER-dist','Convolutive'};
colors={'r','b',[0,180,0]/255};
asz=14;
lw=1.5;
descr=' seizure simulation';

tsc=60; dt=0.05; fs=tsc/dt;
[a b] = butter(3, 2*40/fs, 'low');
[c d] = butter(3,  2*1/fs, 'high');

cols=4;
figure(43); hold on;
sgtitle(['Thin volume',descr],'FontSize',asz+2,'FontWeight','bold')
tstart=5; tend=tstart+10; thigh=tstart+4; tlow=tend-3;

sz_dsS=[-1,-1,-1];
nr_dsS=[-1,-1,-1];

for i=1:3
data=dbase{i}
subplot(4,cols,cols*(i-1)+[1 4]); hold on;
Xf = filter(a, b, data{1});
Xf = filter(c, d, Xf);
ts=data{8};

try
    figure(45); hold on;
    subplot(3,1,i); hold on;
    plot(ts,Xf);
    tinit=1;
    tbegin=tstart-1;
    id_enr=find(ts>tbegin,1,'first');
    id_sz=find(ts>tstart,1,'first');
    id_end=find(ts>tend,1,'first');
    id_init=find(ts>tinit,1,'first');
    
    [sz_ISI,sz_mISI,sz_mf,sz_mdelta]=fISI(Xf(id_sz:id_end),ts(id_sz:id_end),1,1,[1,0]);
    sz_dsS(i)=sz_mdelta;
    [nr_ISI,nr_mISI,nr_mf,nr_mdelta]=fISI(Xf(id_init:id_enr),ts(id_init:id_sz),1,1,[1,0]);
    nr_dsS(i)=nr_mdelta;
    
    
    xlabel('seconds')
    xlim([tinit tbegin]);
    ylim([wdw wup])
    title(['normal simulation (mean amplitude ',num2str(round(nr_mdelta,2)),'\muV, mean ISI ',num2str(round(nr_mISI,3)),' s)'])
    saveas(gcf,[path,'normal_sim.png'])
    
    title(['seizure simulation (mean amplitude ',num2str(round(sz_mdelta,2)),'\muV, mean ISI ',num2str(round(sz_mISI,3)),' s)'])
    xlim([tstart tend]);
    saveas(gcf,[path,'seizure_sim.png'])
catch
    disp('Warning! Failed in computing delta!')
end

figure(43); hold on;
subplot(4,cols,cols*(i-1)+[1 4]); hold on;

plot(ts,Xf,'linewidth',lw,'color',colors{i});
ylim([-50 50]);
xlim([0 16]);
ttl=title(slabel{i});
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0;
ttl.HorizontalAlignment = 'left'; 
set(gca,'FontSize',asz,'FontWeight','bold');
ylabel('\muV','FontWeight','bold');
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
set(gca,'XTick',[]);
end

Tlims=[0,1000];
switcher = @(t) time_step(t/tsc,thigh,tlow)+sment(t/tsc,tstart,thigh,0,1)+sment(t/tsc,tend,tlow,0,1);
Ifun = @(t) time_step(t/tsc,tstart,tend);
tspan=[Tlims(1):dt:Tlims(2)];
Inoise=0.1; Iground=1.0; Iavg=5.0; Ibkg=0.5;
Iplot=((Iavg-Iground)*Ifun(tspan)+Iground).*switcher(tspan)+Ibkg;

subplot(4,cols,4*cols-[3 0]); hold on;
plot(ts,Iplot,'k','linewidth',lw);
ttl=title('Stimulation')
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0;
ttl.HorizontalAlignment = 'left'; 
set(gca,'FontSize',asz,'FontWeight','bold')
xlabel('seconds','FontWeight','bold');
ylabel('nA','FontWeight','bold');
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
xlim([0 16])

function [sz_mdelta,nr_mdelta,ToSave,path]=routine(seed,mode,n,Iavg)
close all
rng(seed)

plt_pks=0;

% %%% Normal state %%%
% a=1; b=3; c=1; d=5; s=4; r = 0.006; x0=-1.6; diff=0.0; I1=1.0-diff; I2=1.0+diff;
%%% Discharge state %%%
%a=1; b=3; c=1; d=5; s=4; r = 0.001; x0=-1.6; diff=0; 
a=1; b=3; c=1; d=5; s=4; r = 0.003; x0=-1.6; diff=0;
%Iavg=2.5; 
%Iavg=1.0;
%tstart=30; tend=70;
%tstart=20; thigh=25; tlow=30; tend=40;
tstart=5; tend=tstart+10; thigh=tstart+4; tlow=tend-3;
%Inoise=0.1; Iground=0.45; Iavg=5.0;
Inoise=0.1; Iground=1.0; 
%Iavg=5.0; 
Ibkg=0.5;
%Inoise=0.1; Iground=0.0; Iavg=5.0;
%slot='0';
slot='0s';
%slot='1';
%slot='2';
%slot='3';
%slot='4';
%mode='Slice6'
%mode='Prob_Dist';
%mode='ER';

bnd_flag='0'; 
%K11=0.5; n=0.05; K21=K11; 
K11=150; 
%n=0.05; 
K21=K11; 
K12=K11*n;
dt=0.05; Tlims=[0,1000];
W0=[0.5;-5;2.45;0;0;0];

% if strcmp(mode,'Celegans')
%     [A,pos,N]=Cel_net(seed);
%     pemp=mean(A(:))
% end

if strcmp(mode,'Slice6')
    [A,pos,N]=Slice_net(seed);
    pemp=mean(A(:))
end

if strcmp(mode,'Slice6D')
    [A,pos,N]=SliceD_net(seed);
    pemp=mean(A(:))
end

if strcmp(mode,'Prob_Dist')
    apar=0.7; bpar=10/400;
    %apar=0.4; bpar=8/2000;
    %apar=0.5; bpar=14/2000;
    %apar=1; bpar=20/2000;
    %apar=0.026; bpar=0;
    [A,pos,N]=Pd_net(seed,apar,bpar);
    pemp=mean(A(:))
end

if strcmp(mode,'ER')
    %apar=0.2; bpar=8/2000;
    %apar=0.5; bpar=14/2000;
    %apar=1; bpar=20/2000;
    apar=0.0438; bpar=0;
    [A,pos,N]=Pd_net(seed,apar,bpar);
    pemp=mean(A(:))
end

NE=round(N*0.9); NI=N-NE;

NW0=[ones(NE,1)*W0(1);ones(NE,1)*W0(2);ones(NE,1)*W0(3);ones(NI,1)*W0(4);ones(NI,1)*W0(5);ones(NI,1)*W0(6)]+1*randn(3*N,1);

figure(1); hold on;
muI=[400 300]/2;
%muI=[250 250 500];
%sigmaIxy=300; 
%sigmaIxy=75; 
sigmaIxy=120; 

path=strcat('MDB_Fig4_Sl/sims_',num2str(Iavg),'/',mode,'/N_',num2str(N),'_n_',num2str(n),'_seed_',num2str(seed),'_I_',num2str(Iavg),'_GND_',num2str(Iground),'_Spos_',num2str(muI(1)),'_',num2str(muI(2)),'_radxy_',num2str(sigmaIxy),'_slot_',slot,bnd_flag,'/')
if ~exist(path, 'dir')
       mkdir(path)
end


M=exp(-((pos(:,1)-muI(1)).^2+(pos(:,2)-muI(2)).^2)/(2*sigmaIxy^2));
I1=(Iavg-Iground).*M(1:NE,1)+Inoise*randn(NE,1); 
I2=(Iavg-Iground).*M((NE+1):N,1)+Inoise*randn(NI,1);

figure(21); hold on;
xp = 0:10:400;
yp = 0:10:300;
[Xp,Yp] = meshgrid(xp,yp);
Mp=(Iavg-Iground).*exp(-((Xp-muI(1)).^2+(Yp-muI(2)).^2)/(2*sigmaIxy^2))+Iground;
contourf(Xp,Yp,Mp);
xlabel('\mum')
ylabel('\mum')
title('Intensity localization')
cplot=colorbar;
cplot.Label.String = 'Intensity';
axis equal
saveas(gcf,[path,'loc.png'])

tic
eps=10^(-4);
%N11=max(sum(A(1:NE,1:NE),2),eps);
%N12=max(sum(A(1:NE,NE+1:N),2),eps);
%N21=max(sum(A(NE+1:N,1:NE),2),eps);

N11=NE; N12=NI; N21=NE;
%Ifun = @(t) time_step(t/20,tstart,tend);
tsc=60;

if slot=='0'
    switcher = @(t) time_step(t/tsc,tstart,tend);
    Ifun = @(t) time_step(t/tsc,tstart,tend);
end
if slot=='0s'
    switcher = @(t) time_step(t/tsc,thigh,tlow)+sment(t/tsc,tstart,thigh,0,1)+sment(t/tsc,tend,tlow,0,1);
    Ifun = @(t) time_step(t/tsc,tstart,tend);
end

if slot=='1'
    switcher = @(t) time_step(t/tsc,tstart,tend);
    Ifun = @(t) time_ramp(t/tsc,tstart,thigh,tlow,tend);
end
if slot=='2'
    period=35;
    switcher = @(t) time_step(t/tsc,tstart,tend+period);
    Ifun = @(t) (time_ramp(t/tsc,tstart,thigh,tlow,tend)+time_ramp(t/tsc,tstart+period,thigh+period,tlow+period,tend+period));
end
if slot=='3'
    period=35;
    switcher = @(t) (time_step(t/tsc,tstart,tend)+time_step(t/tsc,tstart+period,tend+period));
    Ifun = @(t) (time_ramp(t/tsc,tstart,thigh,tlow,tend)+time_ramp(t/tsc,tstart+period,thigh+period,tlow+period,tend+period));
end
if slot=='4'
    fss=3;
    switcher = @(t) slot_4(t,tstart,tend,fss);
    Ifun = @(t) slot_4(t,tstart,tend,fss);
end

tspan=[Tlims(1):dt:Tlims(2)];
dtn=dt*100; tsn=[Tlims(1):dtn:Tlims(2)];
ltn=length(tsn); lt=length(tspan); 
J=1; 
common_noise=Ibkg+0.75*randn(J,ltn);
%Noise=max(0,(1.0*(1+0.1*randn(N,J))/J)*common_noise);
Noise=max(0,ones(N,1)*common_noise);
U=zeros(N,lt);
for i=1:N
    U(i,:)=interp1(tsn,Noise(i,:),tspan,'linear');
end
%U=max(0,1+0.25*randn(N,lt));
Iplot=((Iavg-Iground)*Ifun(tspan)+Iground).*switcher(tspan)+Ibkg;
figure(20); hold on;
plot(tspan/tsc,Iplot);
plot(tspan/tsc,Iground*switcher(tspan)+Ibkg,'--g');
xlabel('time (s)')
ylabel('Intensity')
title('Intensity in time');
drawnow
disp('simulating...')
hw = waitbar(0,'Please wait...');
[ts,W] = ode45(@(t,y) HRnetOde(t,y,Tlims(2),dt,A,NE,NI,U,a,b,c,d,r,s,x0,K11,K12,K21,N11,N12,N21,I1,I2,Ifun,Iground,switcher,hw), tspan, NW0);
close(hw)
disp('processing simulation...')
W=W';
ts=ts/tsc;
time_elapsed=toc

%save('backup.mat')
%load('backup.mat')

vsf=0.2*100;
wup=40; wdw=-40;
X1=mean(W(1:1:NE,:),1)*vsf; 
Y1=mean(W((NE+1):1:2*NE,:),1); 
Z1=mean(W((2*NE+1):1:3*NE,:),1); 

X2=mean(W(3*NE+(1:1:NI),:),1)*vsf; 
Y2=mean(W(3*NE+((NI+1):1:2*NI),:),1); 
Z2=mean(W(3*NE+((2*NI+1):1:3*NI),:),1); 

X=mean(W([1:1:NE,3*NE+(1:1:NI)],:),1)*vsf;

ToSave={X,X1,Y1,Z1,X2,Y2,Z2,ts,mode,Iavg};

%save('backup.mat')
%load('backup.mat')


fs=tsc/dt;
fpass=[0.5 70];

figure(5); hold on;
sgtitle(['Average signals (I=',num2str(Iavg),')'])
subplot(2,1,1); hold on;
title('x')
plot(ts,X);
subplot(2,1,2); hold on;
title('filtered x')
order=3;
%[b,a]=butter(order,fpass/(fs/2), 'bandpass');
[a b] = butter(3, 2*40/fs, 'low');
[c d] = butter(3,  2*1/fs, 'high');
if bnd_flag=='1'
    bmin=-45; bmax=bmin-bmin*256/144
    %Xf=bnd(filter(b,a,X),bmin,bmax);
    Xfr = filter(a, b, X);
    Xfr = filter(c, d, Xfr);
    Xf=bnd(Xfr,bmin,bmax);
    %Xfr=filter(b,a,X);
    %Xf=bnd(bandpass(X,fpass,fs),bmin,bmax);
    %Xfr=bandpass(X,fpass,fs);
    plot(ts,Xfr);
    xlabel('seconds')
    saveas(gcf,[path,'avg.png'])
else
    %Xf=filter(b,a,X);
    Xf = filter(a, b, X);
    Xf = filter(c, d, Xf);
    %Xf=bandpass(X,fpass,fs);
    plot(ts,Xf);
    xlabel('seconds')
    saveas(gcf,[path,'avg.png'])
end



figure(23); hold on;
plot(ts,Xf);
tinit=1;
tbegin=tstart-1;
id_enr=find(ts>tbegin,1,'first');
id_sz=find(ts>tstart,1,'first');
id_end=find(ts>tend,1,'first');
id_init=find(ts>tinit,1,'first');
[sz_ISI,sz_mISI,sz_mf,sz_mdelta]=fISI(Xf(id_sz:id_end),ts(id_sz:id_end),1,1,[1,0]);
sz_mdelta
sz_mISI
[nr_ISI,nr_mISI,nr_mf,nr_mdelta]=fISI(Xf(id_init:id_enr),ts(id_init:id_sz),1,1,[1,0]);
nr_mdelta
nr_mISI

xlabel('seconds')
xlim([tinit tbegin]);
ylim([wdw wup])
title(['normal simulation (mean amplitude ',num2str(round(nr_mdelta,2)),'\muV, mean ISI ',num2str(round(nr_mISI,3)),' s)'])
saveas(gcf,[path,'normal_sim.png'])

title(['seizure simulation (mean amplitude ',num2str(round(sz_mdelta,2)),'\muV, mean ISI ',num2str(round(sz_mISI,3)),' s)'])
xlim([tstart tend]);
saveas(gcf,[path,'seizure_sim.png'])

end

function s=sment(t,tin,tout,S0,S1)
s=time_step(t,min(tin,tout),max(tin,tout)).*(S0+S1*min(1,(t-tin).^2./(tout-tin).^2));
end

function V=slot_4(t,tstart,tend,fsw,tsc)
dT=tend-tstart;
fst=dT*fsw;
NS=fst*2-1;
V=0;
for kf=0:2:(NS-1)
    V=V+time_step(t/tsc,tstart+kf*dT/NS,tstart+(kf+1)*dT/NS);
end
end

function Vpass=bnd(V,bmin,bmax)
Vpass=min(max(V,bmin),bmax);
end

function I=time_step(t,tstart,tend)
I=(t>=tstart).*(t<tend);
end

function I=time_ramp(t,tstart,thigh,tlow,tend)
I=(t>=thigh).*(t<tlow)+(t>=tstart).*(t<thigh).*(t-tstart)./(thigh-tstart)+(t>=tlow).*(t<tend).*(1-(t-tlow)./(tend-tlow));
end

function [A,pos,N]=Pd_net(seed,a,b)
rng(seed)
d=[0:10:2000];
P=pd(d,a,b);

N=550;

figure(14); hold on;
title('probability')
plot(d,P)
xlabel('\mum')
ylabel('probability')

pos=rand(N,3).*[400 300 10];
D = squareform(pdist(pos));
A=(rand(N,N)<pd(D,a,b));

pnet=mean(A(:));

[gin,gout]=find(A>0);

L=[gout,gin];

figure(4); hold on;
title('connection length distribution')
dv=dist_vecL(L,pos);
xlabel('\mum')

[indeg,outdeg]=degsL(L,N);

nsud=20;
figure(3); hold on;
subplot(1,2,1); hold on; title('indegree');
hi=histogram(indeg,nsud,'Normalization','pdf');
% hid=histogram(din,nsud,'Normalization','pdf');
% legend({'sim','data'},'Location','northeast');

subplot(1,2,2); hold on; title('outdegree');
ho=histogram(outdeg,nsud,'Normalization','pdf');
% hout=histogram(dout,nsud,'Normalization','pdf');
% legend({'sim','data'},'Location','northeast');
end

function P=pd(d,a,b)
P=a.*exp(-b.*d);
end

function [A,pos,N]=Slice_net(seed)
rng(seed)
load('GABA.mat');

sl=6;
C=(connMat{sl}>0);
idx=or(sum(C)>0,sum(C')>0);
G=C(idx,idx);
indeg6=sum(G);
outdeg6=sum(G');

figure(1); hold on;
subplot(1,2,1); hold on; title('Indegree distribution');
hi=histogram(indeg6,'Normalization','pdf');

figure(1); hold on;
subplot(1,2,2); hold on; title('Outdegree distribution');
ho=histogram(outdeg6,'Normalization','pdf');

% Nr=max(length(indeg6),length(outdeg6));
delta=0.25; Ek=0.1; noise_factor=3; L=1; pup=1; pdw=0;
Nr=550; 
%delta=0.5; Ek=0.1; noise_factor=3; L=1; pup=1; pdw=0;
%delta=1.5; Ek=1; noise_factor=3; L=1; pup=1; pdw=0;
[Gt,post]=spatial_baSl(Nr,delta,indeg6,outdeg6,Ek,noise_factor,L,pup,pdw);

figure(3); hold on;
sgtitle('Slice fitting')

N=max(Gt(:));
Asim=zeros(N,N);
for i=1:size(Gt,1)
 Asim(Gt(i,2),Gt(i,1))=1;
end

rp=randperm(N,N);
A=Asim(rp,rp);
pos=post(rp,:);
pemp=mean(A(:))

figure(1); hold on;
title('graph plot');
NW=sum(A,2);
sz=round(10+NW/2)*0.5;
scatter3(pos(:,1),pos(:,2),pos(:,3),sz,zeros(N,1),'filled')

axis equal
end

function [A,pos,N]=SliceD_net(seed)
disp(seed)
rng(seed)
load('GABA.mat');

sl=6;
C=(connMat{sl}>0);
idx=or(sum(C)>0,sum(C')>0);
G=C(idx,idx);
indeg6=sum(G);
outdeg6=sum(G');

figure(1); hold on;
subplot(1,2,1); hold on; title('Indegree distribution');
hi=histogram(indeg6,'Normalization','pdf');

figure(1); hold on;
subplot(1,2,2); hold on; title('Outdegree distribution');
ho=histogram(outdeg6,'Normalization','pdf');

% Nr=max(length(indeg6),length(outdeg6));
Nr=550; 
%delta=0.5; Ek=0.1; noise_factor=3; L=1; pup=1; pdw=0;
%delta=1.5; Ek=1; noise_factor=3; L=1; pup=1; pdw=0;
delta=0.5; Ek=0.1; noise_factor=3; bpar=10/400;
[Gt,post]=spatial_baSlD(Nr,delta,indeg6,outdeg6,Ek,noise_factor,bpar);

figure(3); hold on;
sgtitle('Slice fitting')

N=size(post,1);
Asim=zeros(N,N);
for i=1:size(Gt,1)
    if and(Gt(i,2)<=N,Gt(i,1)<=N)
        Asim(Gt(i,2),Gt(i,1))=1;
    end
end

rp=randperm(N,N);
A=Asim(rp,rp);
pos=post(rp,:);
pemp=mean(A(:))

figure(1); hold on;
title('graph plot');
NW=sum(A,2);
sz=round(10+NW/2)*0.5;
gplot3(A,pos,'--','color',[0 0 0 0.1],'linewidth',0.01);
scatter3(pos(:,1),pos(:,2),pos(:,3),sz,zeros(N,1),'filled')

axis equal
end

function dW=HRnetOde(t,W,T,dt,A,NE,NI,U,a,b,c,d,r,s,x0,K11,K12,K21,N11,N12,N21,I1,I2,Ifun,Iground,switcher,hw)

x1=W(1:1:NE); y1=W((NE+1):1:2*NE); z1=W((2*NE+1):1:3*NE); 
x2=W(3*NE+(1:1:NI)); y2=W(3*NE+((NI+1):1:2*NI)); z2=W(3*NE+((2*NI+1):1:3*NI)); 
%X1=mean(x1); X2=mean(x2);
N=NE+NI; 

dx1=y1-a.*x1.^3+b*x1.^2-z1+(K11.*sum(A(1:NE,1:NE).*(x1'-x1),2)./N11-K12.*sum(A(1:NE,NE+1:N).*(x2'-x1),2)./N12)+(I1.*Ifun(t)+Iground)*switcher(t)+U(1:NE,round(t/dt)+1);
dy1=c-d.*x1.^2-y1;
dz1=r.*(s.*(x1-x0)-z1);%-0.1*z1.^7.*(z1<0));

dx2=y2-a.*x2.^3+b.*x2.^2-z2+K21.*sum(A(NE+1:N,1:NE).*(x1'-x2),2)./N21+(I2.*Ifun(t)+Iground)*switcher(t)+U(NE+1:N,round(t/dt)+1);
dy2=c-d.*x2.^2-y2;
dz2=r.*(s.*(x2-x0)-z2);%-0.1*z2.^7.*(z2<0));

dW=[dx1; dy1; dz1; dx2; dy2; dz2];

waitbar(t/T,hw)

end


