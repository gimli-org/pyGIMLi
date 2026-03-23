function writeroundpoly(Shotpos,filename,dd,iseasy)

% WRITEPOLY - Creatte poly file from Shot positions
% writepoly(Shotpos,filename[,refine,easymesh])
% Shotpos ... Nx2 matrix of positions
% filename .. File name for poly output
% refine  ... relative refinement near positions
% iseasy  ... write easy mesh instead of triangle

if nargin<1, error('No shot provided'); end
if nargin<2, filename='test.poly'; end
if nargin<3, dd=0.3; end
if nargin<4, iseasy=0; end
if size(Shotpos,2)<2, Shotpos(:,2)=0; end

np=size(Shotpos,1);
mid=mean(Shotpos);
relpos=Shotpos-repmat(mean(Shotpos),np,1);
rad=sqrt(sum(relpos.^2,2));
ang=atan2(relpos(:,2),relpos(:,1));
sang=sort(ang);
dpos=1;
% if iseasy||(dd<=0)||(dd>=1), % no refinement
if (dd<=0)||(dd>=1), % no refinement
    newang=sang;
elseif dd>0.33, %treat as 0.5 (add one more)
    dang=mod([diff(sang);sang(1)-sang(end)],2*pi);
    newang=reshape([sang sang+dang/2]',np*2,1);    
    dpos=2;
else % real refinement
    dang=mod([diff(sang);sang(1)-sang(end)],2*pi);
    newang=reshape([sang-dang*dd sang sang+dang*dd]',np*3,1);
    dpos=3;
end
newrad=interp1(ang,rad,newang,'linear','extrap');
%spez. Behandlung von newrad(1) und newrad(end)?
newpos=[newrad.*cos(newang) newrad.*sin(newang)]+repmat(mean(Shotpos),length(newang),1);
np=size(newpos,1);
% nn=[(0:np-1)' newpos ones(np,1)*(-99)];
nn=[(0:np-1)' newpos ones(np,1)*(0)];
% nn(1+1:dpos:end,4) = -99;
nn(1+(dpos>2):dpos:end,4) = -99;
ee=nn;
ee(:,2)=ee(:,1);ee(:,4)=-1;
ee(:,3)=ee(:,2)+1;ee(end,3)=0;
if iseasy,
    fid=fopen(filename,'w');
    fprintf(fid,'%d\n',np);    
    nn(:,4)=dd;
    nn(:,5)=1;
    fprintf(fid,'%d:\t%f\t%f\t%f\t%d\n',nn');
    fprintf(fid,'%d\n',np);
    ee(:,4)=1;
    fprintf(fid,'%d:\t%d\t%d\t%d\n',ee');
    fclose(fid);
else
    fid=fopen(filename,'w');
    fprintf(fid,'%d %d %d %d\n',np,2,0,1);
    fprintf(fid,'%d\t%.3f\t%.3f\t%d\n',nn');
    fprintf(fid,'%d %d\n',np,1);
    fprintf(fid,'%d\t%d\t%d\t%d\n',ee');
    fprintf(fid,'%d\n%d\n',0,0);
    fclose(fid);
end