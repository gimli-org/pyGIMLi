function [cols,cmin,cmax,cmap,islog]=getcols(att,MAL,alfa)

% GETCOLS - Returns Color (IxJx3) matrix for patch plots
% cols=getcols(field,spec,alfa)

if (nargin<3)||(~isequal(numel(alfa),numel(att))), 
    alfa=ones(size(att)); end
if nargin<2, MAL=[]; end
cols=zeros([size(att) 3]);

% if ~isfield(MAL,'cauto')&isfield(MAL,'cmin')&isfield(MAL,'cmax'), MAL.cauto=0; end
if ~isfield(MAL,'cauto'), MAL.cauto=~(isfield(MAL,'cmin')&isfield(MAL,'cmax')); end
if isfield(MAL,'cmin'), cmin=MAL.cmin; else cmin=min(att); end
if isfield(MAL,'cmax'), cmax=MAL.cmax; else cmax=max(att); end
if isfield(MAL,'clog'), islog=MAL.clog;
else islog=(min(att)>0)&&(cmin>0); end
if islog, 
    att=log10(att); 
    cmin=log10(cmin);
    cmax=log10(cmax);
end
perc=5;if isfield(MAL,'perc'), perc=MAL.perc; end
if (~isfield(MAL,'cauto'))||(MAL.cauto>0), 
    if (length(unique(att))<10), cmin=min(att);cmax=max(att); 
    else [cmin,cmax]=interperc(att,[perc 100-perc]); end
end
cmap=jet(64); 
if (~isfield(MAL,'cmap')||(MAL.cmap~=0))&&((~islog)&&(max(att)*min(att)<0)),
    cmap=colormap(b2r(64));
    if ~isfield(MAL,'cauto')||(MAL.cauto==1), 
        cmax=max(abs([cmin cmax]));cmin=-cmax; end
else
    if isfield(MAL,'cmap'),
        if MAL.cmap==2, cmap=b2r; end
        if MAL.cmap==3, cmap=gray(64); end
        if MAL.cmap==5, hh=hsv(6);cmap=hh(5:-1:1,:); end
%         if MAL.cmap==6, cmap=hsv(6); end
        if MAL.cmap==6, hh=hsv(6);cmap=hh([5:-1:1 6],:); end
        if MAL.cmap==7, cmap=[0 0 1;0 1 1;0 1 0;1 1 0;1 0 0;1 0 1;0.5 0.5 0.5]; end
    end
end
if isfield(MAL,'cflip')&&(MAL.cflip), 
    cmap=flipud(cmap);colormap(cmap);
end
if ~(cmax>cmin), cmax=cmin+1; end
lcm=length(cmap);
cind=round(1+(att-cmin)/(cmax-cmin)*(lcm-1));
cind(cind<1)=1;
cind(cind>lcm)=lcm;
s1=size(att,1);
for k=1:3,
   cols(:,:,k)=cmap(cind(:),k).*alfa(:)+1-alfa(:);
end
return
for i=1:size(att,1),
    for j=1:size(att,2),
        cols(i,j,:)=cmap(cind(i,j),:)*alfa(i,j)+1-alfa(i,j);
    end
end