function [L,Cx,Cz] = smooth2d1st(x,z,nor,az)

% SMOOTHMAT2D - Computes smoothness matrix
% L = smooth2d1st(x,z)
% L = smooth2d1st(x,z,nor,az)
% [L,Cx,Cz] = smooth2d1st(...)
% L constructed by Cx'*Cx+Cz'*Cz
% or, equivalently, L=C'*C with C=[Cx Cz]
% x,z .. model grid lines
% nor .. normalize vectors x and y (default=0)
% az  .. factor for weighting z-derivarives (default=1)

if nargin<4, az=1; end
if nargin<3, nor=0; end
argmin=1;
if min(size(x))>1, % Matrix
    z=0:size(x,2);
    x=0:size(x,1);
else
    argmin=2;
end
if nargin<argmin, error('Too less input arguments!'); end
if nor, x=1:length(x);z=1:length(z); end
nx=length(x)-1;nz=length(z)-1;
x=(x(1:end-1)+x(2:end))/2; % Midpoint
z=(z(1:end-1)+z(2:end))/2;
xmed=median(diff(x)); % mean diff(x)
x=x/xmed;z=z/xmed;
nn=nx*nz;
one=repmat([1./diff(x(:));0],nz,1);
Cx=spdiags([-one one],[0 1],nn-1,nn);
Cx(nx:nx:end,:)=[];
two=repmat(1./diff(z(:))',nx,1);%two(:,3)=0; %
Cz=spdiags([-two(:) two(:)],[0 nx],nx*(nz-1),nn);
Cz=sqrt(az)*Cz;
L=Cx'*Cx+Cz'*Cz;