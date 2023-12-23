function [ISI,mISI,mf,mdelta]=fISI(Xs,t,cutoff,DV,plt)

[npks,nidxs] = findpeaks(-Xs); npks=-npks;
[ppks,pidxs] = findpeaks(Xs);
lp=min(length(nidxs), length(pidxs));
ISI=abs(t(nidxs(1:lp))-t(pidxs(1:lp)));
delta=abs(npks(1:lp)-ppks(1:lp));
valid=find(and(ISI(:)<cutoff,delta(:)>DV));
ISI=ISI(valid);
mISI=mean(ISI);
mdelta=mean(delta(valid));
mf=1/mISI;

if plt(1)==1
    plot(t(nidxs(valid)),npks(valid),'rx')
    plot(t(pidxs(valid)),ppks(valid),'gx')
end

if plt(2)==1
    figure(2); hold on;
    histogram(ISI)
end
end