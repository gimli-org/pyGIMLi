function Poly = create2dpoly(Shotpos,zz,dd,zvec)

% WRITEPOLY - Creatte poly file from Shot positions
% writepoly(Shotpos,filename[,maxdepth,refinement])

if nargin<1, error('No shot provided'); end
if nargin<2, zz=(max(Shotpos(:,1))-min(Shotpos(:,1)))/3; end
if nargin<3, dd=0; end
if nargin<4, zvec=[]; end
if size(Shotpos,2)<2, Shotpos(:,2)=0; end

dx=median(diff(Shotpos(:,1)));
di=1;
if dd>0, % refinement     
    newpos=Shotpos;
    if dd>0.3, % set to 0.5
        newpos(1:end-1,2)=newpos(1:end-1,1)+diff(newpos(:,1))/2;
        newpos=reshape(newpos',[],1);newpos(end)=[];
        di=2;
    else % right and left hand side
        newpos(:,2)=newpos(:,1);
        dnew=diff(newpos(:,1))*dd;
        newpos(2:end,1)=newpos(2:end,2)-dnew;
        newpos(1:end-1,3)=newpos(1:end-1,2)+dnew;
        newpos=reshape(newpos',[],1);
        newpos([1 end])=[];
        di=3;
    end
    newpos(:,2)=interp1(Shotpos(:,1),Shotpos(:,2),newpos,'spline');
    Shotpos=newpos;
end
mi=min(Shotpos(:,1))-dx;
ma=max(Shotpos(:,1))+dx;
if nargin<3, zz=(ma-mi)/4; end
np=size(Shotpos,1);
nn=(0:np+3)';
nn(:,2)=[mi;Shotpos(:,1);ma;ma;mi];
nn(:,3)=[Shotpos(1,2);Shotpos(:,2);Shotpos(end,2);Shotpos(end,2);Shotpos(1,2)];
nn(end-1:end,3)=nn(end-1:end,3)-abs(zz);
nn(2:di:end-2,4)=-99;
nn([1 end-2:end],4)=0;
ee=nn;
ee(:,2)=ee(:,1);ee(:,4)=-1;
ee(:,3)=ee(:,2)+1;ee(end,3)=0;
if length(zvec)>0,
   ln=size(nn,1);
   le=size(ee,1);
   lz=length(zvec);
   apos=nn(1,2:3); % first surface point
   epos=nn(ln-2,2:3); % last surface point
   for i=1:length(zvec), %add new points
       nn(ln+i,2:3)=epos-[0 1]*zvec(i); % right
       nn(ln+lz+i,2:3)=apos-[0 1]*zvec(i); % left
   end
   nn(:,1)=0:size(nn,1)-1;
   ee(end-2:2:end,:)=[]; % delete old vertical lines
   re=ln-3;li=0;
   for i=1:length(zvec),
      ee(le+i-2,2:3)=[re ln-1+i]; 
      ee(le+lz+i-2,2:3)=[ln+lz-1+i li];
      li=ln+lz-1+i;re=ln-1+i;
   end
   ee(end+1,2:3)=[re ln-2];
   ee(end+1,2:3)=[ln-1 li];
   for i=1:length(zvec),
      ee(end+1,2:3)=[ln-1+i ln+lz-1+i]; 
   end
   ee(:,1)=0:size(ee,1)-1;
   ee(:,4)=-1;
end
Poly.node=nn(:,2:end);
Poly.edge=ee(:,2:end)+1;