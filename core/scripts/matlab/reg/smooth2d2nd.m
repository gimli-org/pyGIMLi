function C=smooth2d2nd(x,z,rb,nor,az)

% SMOOTH2D2nd - Second order smoothness matrix for 2d models
% C = smooth2d2nd(x,z[,bc[,nor[,alfaz]]])
% x,z   - vectors of grid lines
% bc    - boundary condition index (default=0)
%        (0=dirichlet, 1=neumann, 2=modified neumann)
% nor   - normalize vectors (default=0)
% alfaz - weighting factor for z-direction

if nargin<2, % for testing
    x=[0 1 2 4];x=0:3;
    z=x; 
end
if nargin<3, rb=0; end
if nargin<4, nor=0; end
if nargin<5, az=1; end
if nor,
    x=1:length(x);z=1:length(z);
end
% midpoints representing model cells
xq=(x(1:end-1)+x(2:end))/2;
zq=(z(1:end-1)+z(2:end))/2;
% scaling because of comparability of regularization parameter
xmed=median(diff(xq));
xq=xq/xmed;
zq=zq/xmed;
dx=diff(xq(:));
dx=[dx(1);dx;dx(end)];
lx=length(dx)-1;
dz=diff(zq(:));
dz=[dz(1);dz;dz(end)];
lz=length(dz)-1;
[DX,DZ]=ndgrid(dx,dz);
DX(:,end)=[];
DZ(end,:)=[];
DXQ=DX(1:end-1,:)+DX(2:end,:);
DZQ=DZ(:,1:end-1)+DZ(:,2:end);
im=1./DX(1:end-1,:)./DXQ;
ip=1./DX(2:end,:)./DXQ;
km=az./DZ(:,1:end-1)./DZQ;
kp=az./DZ(:,2:end)./DZQ;
if rb==0, ik=-im-ip-km-kp; end % dirichlet
im(1,:)=0;
ip(end,:)=0;
km(:,1)=0;
kp(:,end)=0;
if rb>1, % mod neumann
    im(end,:)=0;
    ip(1,:)=0;
    km(:,end)=0;
    kp(:,1)=0;
end
if rb>0, ik=-im-ip-km-kp; end
lxz=lx*lz;
C=spdiags([km(:) im(:) ik(:) ip(:) kp(:)],...
    [lx 1 0 -1 -lx],lxz,lxz)';