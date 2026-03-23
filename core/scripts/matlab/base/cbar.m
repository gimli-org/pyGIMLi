function cb=cbar(cmin,cmax,lolo,dir,anz,loesch,anot)

% CBAR - Draw single color bar (for exporting)
% cbar(cmin,cmax[,logspaced,dir,nticks,delfirst,annotation])
% cmin/cmax - minimum/maximum color
% logspaced - logarithmic spacing
% dir - direction (0=horizontal,1=vertical)
% nticks - number of ticks
% delfirst - delete first (for alpha blending)
% anotation - annotation (axis unit)

if (nargin==1)&&isstruct(cmin),
    lolo=0;anz=7;loesch=0;dir=0;
    mal=cmin;
    cmin=mal.cmin;cmax=mal.cmax;
    if isfield(mal,'clog'), lolo=mal.clog; end
    if isfield(mal,'canot'), anot=mal.canot; end
else
    if nargin<7, anot=''; end
    if nargin<6, loesch=0; end
    if nargin<5, anz=0; end
    if nargin<4, dir=0; end %horiz
    if nargin<3, lolo=0; end
    if nargin<2, cmax=cmin+1; end
    if nargin<1, cmin=0; end
end
if lolo, cmax=log10(cmax);cmin=log10(cmin); end
n=size(colormap,1);
if anz==0, anz=11-dir*4; end
if dir==0,
    image([1 n 1],[0 n/10],(1:n));
    cb=gca;
%     axis tight
    set(cb,'XLim',[0.5 n+0.5]);
    set(cb,'DataAspectRatio',[1 10 1]);
    set(gca,'YTick',[],'XDir','normal');
    yt=linspace(loesch,n,anz);
    set(gca,'XTick',yt+0.5);
    ytl=yt/n*(cmax-cmin)+cmin;
    if lolo==1, ytl=10.^ytl; end
    fi=find(ytl>=9.5);ytl(fi)=round(ytl(fi));
    fi=find(ytl<9.5);
    for i=1:length(fi),
        l=0;yy=ytl(fi(i));
        if yy~=0,
            while abs(yy)<9.5, yy=yy*10;l=l+1; end
            yy=round(yy);
            for ll=1:l, yy=yy/10; end
        end
        ytl(fi(i))=yy;
    end
    set(gca,'XTickLabel',ytl);
    if ~isempty(anot),
        xl=get(gca,'Xlim');yl=get(gca,'Ylim');
%         set(text(mean(xl),yl(1),anot),'VerticalAlignment','bottom','HorizontalAlignment','center');
        set(text(xl(end),mean(yl),[' ' anot]),'VerticalAlignment','middle','HorizontalAlignment','left');
    end
    %     po=get(gca,'Position');
%     po(4)=po(4)/2;
%     po(2)=po(2)*2;
%     set(gca,'Position',po);
else
    image([0 n/10],[1 n 1],(1:n)');
    set(gca,'DataAspectRatio',[5 1 1]);
%     axis tight
    set(gca,'YLim',[0.5 n+0.5]);
    set(gca,'XTick',[],'YDir','normal');
    yt=linspace(loesch,n,anz);
    set(gca,'YTick',yt+0.5);
    ytl=yt/n*(cmax-cmin)+cmin;
    if lolo==1, ytl=10.^ytl; end
    fi=find(ytl>=9.5);ytl(fi)=round(ytl(fi));
    fi=find(ytl<9.5);
    for i=1:length(fi),
        l=0;yy=ytl(fi(i));
        if yy~=0,
            while abs(yy)<9.5, yy=yy*10;l=l+1; end
            yy=round(yy);
            for ll=1:l, yy=yy/10; end
        end
        ytl(fi(i))=yy;
    end
    set(gca,'YTickLabel',ytl,'YAxisLocation','right');
    if ~isempty(anot),
        xl=get(gca,'Xlim');yl=get(gca,'Ylim');
%         t=text(xl(end),mean(yl),anot);
%         set(t,'VerticalAlignment','middle','HorizontalAlignment','right');
        t=text(mean(xl),yl(end)+1,anot);
        set(t,'VerticalAlignment','bottom','HorizontalAlignment','center');
    end
end
% if cmin>0,
%     title('\rho in \Omega\cdotm');
% end