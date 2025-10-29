function OUT=minmaxfwd(res,thk,coilsep,freq,zz)

% MINMAXFWD - Compute FDEM response for separated coils (MINMAX)
% OUT=minmaxfwd(res,thk,coilsep,freq,zz)

if nargin<5, zz=-1; end
if nargin<4, freq=[110,220,440,880,1760,3520,7040,14080,28160,56320]; end
if nargin<3, coilsep=50; end
if nargin<2, thk=10; end
if nargin<1, res=[1 1]; end

%%
ze=-abs(zz);zs=ze; % both at same height
rm=sqrt(coilsep.^2+(ze-zs).^2);
rp=sqrt(coilsep.^2+(ze+zs).^2);
fas=-(3*(ze+zs).^2-rp.^2)./rp.^5/4/pi;
%ff=(3*(ze-zs).^2-rp.^2)./rp.^5; % Hz
%bh=1/(4*pi)+0*100; % A/m oder nT
%fas=-bh*ff;

%%
nf=length(freq);
OUT=zeros(nf,2);
for i=1:nf,
    fields=vmd_f([1 0 0],freq(i),res,thk,-abs(zz),-abs(zz),coilsep,1,1);
    OUT(i,1)=real(fields(1));
    OUT(i,2)=imag(fields(1));
end
% if length(fas)==length(freq),
   OUT(:,1)=-OUT(:,1)./fas(:)*100-100;
   OUT(:,2)=-OUT(:,2)./fas(:)*100;
% end