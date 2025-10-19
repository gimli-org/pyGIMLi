function misfit=plotmrsdatafit(t,Data,Resp,cmin,cmax)

if nargin<4, cmin=0; end
if nargin<5, cmax=0; end

if cmin==0, cmin=interperc(log10(abs(Data(:))),3); end
if cmax==0, cmax=interperc(log10(abs(Data(:))),97); end

nT=size(Data,2);
subplot(3,3,[1 4]);
imagesc(t,1:nT,log10(abs(Data')));
caxis([cmin cmax]);colorbar horiz
subplot(3,3,[2 5]);
imagesc(t,1:nT,log10(abs(Resp')))
caxis([cmin cmax]);colorbar horiz
subplot(3,3,[3 6]);
misfit=Data-abs(Resp);
imagesc(t,1:nT,misfit');colorbar horiz
subplot(3,3,7:9);hist(misfit(:),100);

