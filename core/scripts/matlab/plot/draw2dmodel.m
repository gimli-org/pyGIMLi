function [cmin,cmax]=draw2dmodel(x,z,dM,MAL,Alpha)

% DRAW2DMODEL - Draw 2D Model
% draw2dmodel(x,z,Model[,OPT])
% x,z..model block corners
% dM..Model(length(x)-1-2*border,length(z)-1-border)
% OPT..structure of possible fields
%      xdir - direction of  x axis (0/1)
%      cauto - automatic coloring (0/1)
%      cmin/cmax - minimum/maximum color values
%      cmap  - colormap (0-6 default,b2r,hot,gray,jet,cool)
%      log   - logarithmic coloring
%      style - Color boxes(0), filled contours(1), interpolated(2)
%      high  - plot z-axis with factor x
%Alpha=[];

global N
if nargin<4,
    MAL=struct('cauto',1,'cmin',100,'cmax',500,'cmap',0,'log',0,'xdir',0,'style',0,'high',1);
end
if nargin<5, Alpha=1; end
if isempty(Alpha), Alpha=1; end
clog=0;
if isfield(MAL,'xdir'), xdir=MAL.xdir; else xdir=0; end
if isfield(MAL,'cmax'), cmax=MAL.cmax; else cmax=1; end
if isfield(MAL,'cmin'), cmin=MAL.cmin; else cmin=0; end
if isfield(MAL,'cmap'), cmap=MAL.cmap; else cmap=0; end
if isfield(MAL,'clog'), clog=MAL.clog; end
if isfield(MAL,'log'), clog=MAL.log; end
if isfield(MAL,'cauto'), cauto=MAL.cauto; else cauto=1; end
if isfield(MAL,'cbar'), cbar=MAL.cbar; else cbar=1; end
if isfield(MAL,'elec'), elec=MAL.elec; else elec=0; end
if isfield(MAL,'high'), high=MAL.high; else high=1; end
if isfield(MAL,'style'), style=MAL.style; else style=0; end
if isfield(MAL,'cont'), cont=MAL.cont; else cont=15; end
if isfield(MAL,'alpha')&&(MAL.alpha==0), Alpha=0; end
if (cauto==0)&&(clog==1), cmin=log10(cmin);cmax=log10(cmax); end
cont=30;
if max(dM(:))==min(dM(:)), style=0; end
% if min(dM(:))*max(dM(:))<0, % symmetric blue-white-red
%     cmax=max(abs(dM(:)));
%     cmin=-cmax;
%     cmap=2;
%     clog=0;
%     cauto=0;
% end

if clog==1,
    dm=log10(reshape(dM,length(x)-1,length(z)-1))';
else
    dm=reshape(dM,length(x)-1,length(z)-1)';
end
uni=unique(dm(:));
if length(uni)<10, Alpha=1; end
if cauto==1,
    if length(uni)<10,
        cmin=min(uni);
        cmax=max(uni);
    else
        [NN,VV]= hist(dm(:),100);
        CN=cumsum(NN);
        CN=CN/max(CN);
        imin=max(find(CN<0.01));
        imax=min(find(CN>0.99));
        if isempty(imin), imin=1; end
        if isempty(imax), imax=length(VV); end
        cmin=VV(imin);
        cmax=VV(imax);
    end
    if cmax<=cmin, cmax=max(dm(:));cmin=min(dm(:)); end
    if cmax<=cmin, cmax=cmin*1.01;cmin=cmin*0.99; end
end
if isequal(size(Alpha),size(dM)), % Alpha shading
    if isfield(MAL,'noascal')&&(MAL.noascal),
        AA=Alpha';
    else
        AA=log(Alpha');
        [nn,hh]=hist(AA(:),50);
        nnn=cumsum(nn)/prod(size(AA));
        mi=hh(min(find(nnn>0.1)));
        ma=hh(max(find(nnn<0.7)));
        AA=(AA-mi)/(ma-mi);
        AA(AA<0)=0;
        AA(AA>1)=1;
    end
end
if style==1, %(>0) filled contours( or smooth surface)
    xx=(x(1:end-1)+x(2:end))/2;
    zz=(z(1:end-1)+z(2:end))/2;
    zz=[z(1);zz(:);z(end)];
    xx=[x(1);xx(:);x(end)];
    dm=[dm(1,:);dm;dm(end,:)];
    dm=[dm(:,1),dm,dm(:,end)];
    if style>1, %smooth surface
        pcolor(xx,zz,dm);
        shading flat
        grid off
        if size(Alpha)==size(dM),
            AA=[AA(1,:);AA;AA(end,:)];
            AA=[AA(:,1:end),AA(:,end),AA(:,end)];
            alpha(AA);
        end
        shading interp;
    else % filled contours
        contourf(xx,zz,dm,cont);
        shading flat
    end
else % colorboard
    if size(Alpha)==size(dM), % Alpha shading
        vv=version;
        if str2double(vv(1:3))<6.5, % altes matlab
            dm(end+1,end+2)=10000;
            pcolor([x(1);x(:)],z,dm);
            shading flat
            AA(end,end+1)=0;
            alpha(AA);
        else % matlab version 6.5(R13)
            dm(end+1,:)=dm(end,:);
            dm=[dm,dm(:,end)];
            pcolor(x,z,dm);
            AA(end+1,end+1)=1;
            alpha(AA(:,2:end));
        end
    else %alpha
        dm(end+1,:)=dm(end,:);
        dm(:,end+1)=dm(:,end);
        pcolor(x(1:end),z(1:end),dm);
    end %alpha
    if (length(x)>30)&&(length(unique(dM))>5),shading flat; end
end

if high>0, 
%     axis equal; 
set(gca,'DataAspectRatio',[1 1 1]);
end
axis tight
set(gca,'XAxisLocation','top');
set(gca,'Ydir','reverse');
if xdir>0,
    set(gca,'XDir','reverse'); 
else
    set(gca,'XDir','normal');
end
switch cmap
    case 1,
        colormap default
    case 2,
        colormap(b2r);
    case 3,
        colormap(hot);
    case 4,
        colormap(gray);
    case 5,
        colormap(jet);
    case 6,
        colormap(cool);
    otherwise,
        colormap default
end
caxis([cmin,cmax]);
if cbar,
    hc=colorbar('horiz');
    if(clog==1)
        clabel=10.^get(hc,'XTick');
        fi=find(clabel>10);
        clabel(fi)=round(clabel(fi));
        fi=find((clabel<=10)&(clabel>1));
        clabel(fi)=round(10*clabel(fi))/10;
        fi=find((clabel<1)&(clabel>=0.01));
        clabel(fi)=round(100*clabel(fi))/100;
        set(hc,'XTickLabel',num2strcell(clabel));
%         set(hc,'XTickLabel',num2strcell(10.^get(hc,'XTick')));
        %set(hc,'YTick',mean(get(hc,'YLim')));
        %set(hc,'YTickLabel','R in Ohmm');
        %set(hc,'XLabel','R in Ohmm');
        %set(get(hc,'XLabel'),'String','\rho in \Omega m');
    end
end
if high>1,
    set(gca,'DataAspectRatio',[1 1/high 1]);
end
if length(z)>6,
    if size(z,1)==1,
        zz=[z(1) z(3) z(5:end)];
    else
        zz=[z(1);z(3);z(5:end)];
    end;
else
    zz=z;
end
zz=zz(:);
xl=get(gca,'XTickLabel');
if size(xl,1)>1, xl(end-1,1:3)='x/m'; end
set(gca,'XTickLabel',xl);
ddz=round(length(zz)/5);if ddz==0, ddz=1; end
set(gca,'YTick',[zz(1:ddz:end-1);zz(end)]);
yt=get(gca,'YTick');
yt(end:end+1)=[(yt(end-1)+yt(end))/2 yt(end)];
set(gca,'YTick',yt);
% yl=get(gca,'YTickLabel');yl(end-1,1:3)='z/m';yl(end-1,4:end)=' ';
yl=num2strcell(yt);yl{end-1}='z/m';
set(gca,'YTickLabel',yl);

if elec,
    if isfield(N,'elec'),
        hold on
        plot(N.elec(:,1),N.elec(:,2),'kv');
        hold off
    end
else
%     hold on
%     plot(x(1),0,'.','MarkerSize',1);
%     hold off
end
if style==2, shading('interp'); end
lx=length(x);lz=length(z);
hold on;set(line(x([1 lx lx 1 1]),z([1 1 lz lz 1])),'Color','black');hold off;
%Testversion n. 5 Zeilen einklammern
%      set(line(x([1 end]),z([1 end])),'Color','black');
%      set(line(x([1 end]),z([end 1])),'Color','black');
%      tv=[145 144 150 140 141 154 169 223 139 140 154 171];
%      tt=text(mean(x),mean(z),char(255-fliplr(tv)));
%      set(tt,'FontSize',24,'HorizontalAlignment','center','VerticalAlignment','middle');
if clog, cmin=10^cmin;cmax=10^cmax; end
if nargout==1, cmin=[cmin cmax]; end