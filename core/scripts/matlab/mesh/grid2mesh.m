function [Mesh,T] = grid2mesh(x,z,M)

% GRID2MESH - Transform regular (2d) grid to triangle mesh
%             by dividing each rectangle into 2 triangles
% [Mesh,T] = grid2mesh(x,z[,M])

if nargin<2, error('Too less input arguments (x,y)!'); end
if nargin<3, M=ones(length(x)-1,length(z)-1); end

Mesh=[];Mesh.dim=2;
nx=length(x);nz=length(z);
x=x(:);z=z(:);
dx=diff(x);dz=diff(z);
Mesh.nnodes=(nz-1)*(2*nx-1)+nx;
Mesh.node=zeros(Mesh.nnodes,Mesh.dim);
Mesh.node(1:nx,1)=x;Mesh.node(1:nx,2)=z(1);
Mesh.nodemarker=zeros(Mesh.nnodes,1);
Mesh.ncells=4*(nx-1)*(nz-1);
Mesh.cell=zeros(Mesh.ncells,Mesh.dim+1);
Mesh.cellnodes=ones(Mesh.ncells,1)*(Mesh.dim+1);
Mesh.nbounds=6*(nx-1)*(nz-1)+nz-1+nx-1;
Mesh.boundnodes=ones(Mesh.nbounds,1)*2;
Mesh.boundmarker=zeros(Mesh.nbounds,1);
Mesh.bound=zeros(Mesh.nbounds,Mesh.dim);
Mesh.bound(1:nx-1,1)=(1:nx-1)';Mesh.bound(1:nx-1,2)=(2:nx)';
ib=nx-1;in=nx;ic=0; % bound, nodes and cell counter
nm=numel(M);T=spalloc(nm*4,nm,nm*4);iii=[0 nx-1 2*nx-2 3*nx-3];
for iz=1:nz-1, % all layers
    os=(iz-1)*4*(nx-1);
    for i=1:nx-1,
        ii=(iz-1)*(nx-1)+i;
        T(os+i+iii,ii)=1;
    end
    % 2 layers of nodes
    Mesh.node(in+1:in+nx-1,1)=x(1:end-1)+dx/2;
    Mesh.node(in+1:in+nx-1,2)=z(iz)+dz(iz)/2;
    Mesh.node(in+nx:in+2*nx-1,1)=x(:);
    Mesh.node(in+nx:in+2*nx-1,2)=z(iz+1);
    lo=(in-nx+1:in-1)';
    ro=(in-nx+2:in)';
    mm=(in+1:in+nx-1)';
    lun=(in+nx:in+2*nx-2)';
    ru=(in+nx+1:in+2*nx-1)';
    % upper triangle
    Mesh.cell(ic+1:ic+nx-1,1)=lo;
    Mesh.cell(ic+1:ic+nx-1,2)=ro;
    Mesh.cell(ic+1:ic+nx-1,3)=mm;
    Mesh.cellattr(ic+1:ic+nx-1)=M(:,iz);
    ic=ic+nx-1;
    % left triangle
    Mesh.cell(ic+1:ic+nx-1,1)=lo;
    Mesh.cell(ic+1:ic+nx-1,2)=mm;
    Mesh.cell(ic+1:ic+nx-1,3)=lun;
    Mesh.cellattr(ic+1:ic+nx-1)=M(:,iz);
    ic=ic+nx-1;
    % right triangle
    Mesh.cell(ic+1:ic+nx-1,1)=ro;
    Mesh.cell(ic+1:ic+nx-1,2)=ru;
    Mesh.cell(ic+1:ic+nx-1,3)=mm;
    Mesh.cellattr(ic+1:ic+nx-1)=M(:,iz);
    ic=ic+nx-1;
    % lower triangle
    Mesh.cell(ic+1:ic+nx-1,1)=lun;
    Mesh.cell(ic+1:ic+nx-1,2)=mm;
    Mesh.cell(ic+1:ic+nx-1,3)=ru;
    Mesh.cellattr(ic+1:ic+nx-1)=M(:,iz);
    ic=ic+nx-1;
    Mesh.bound(ib+1:ib+nx,1)=(in-nx+1:in)'; %vertical line
    Mesh.bound(ib+1:ib+nx,2)=(in+nx:in+2*nx-1)';
    ib=ib+nx;
    Mesh.bound(ib+1:ib+nx-1,1)=lo;Mesh.bound(ib+1:ib+nx-1,2)=mm;ib=ib+nx-1;
    Mesh.bound(ib+1:ib+nx-1,1)=ro;Mesh.bound(ib+1:ib+nx-1,2)=mm;ib=ib+nx-1;
    Mesh.bound(ib+1:ib+nx-1,1)=lun;Mesh.bound(ib+1:ib+nx-1,2)=mm;ib=ib+nx-1;
    Mesh.bound(ib+1:ib+nx-1,1)=ru;Mesh.bound(ib+1:ib+nx-1,2)=mm;ib=ib+nx-1;
    Mesh.bound(ib+1:ib+nx-1,1)=lun;Mesh.bound(ib+1:ib+nx-1,2)=ru;ib=ib+nx-1;
    in=in+2*nx-1;
end
% aa=(0:30)';aa(:,2)=0;fi=find(ismember(Mesh.node,aa,'Rows'));
% Mesh.nodemarker(fi)=-99;
Mesh.node(:,end)=-Mesh.node(:,end);
Mesh.cellattr=Mesh.cellattr(:);
Mesh.boundleft=zeros(Mesh.nbounds,1);
Mesh.boundright=Mesh.boundleft;
for i=1:Mesh.nbounds,
    bb=Mesh.bound(i,:); % edge nodes
    [fi,jj]=find(Mesh.cell==bb(1));
    for j=2:length(bb),
        [ii,jj]=find(Mesh.cell==bb(j));
        fi=intersect(fi,ii);
    end
    Mesh.boundleft(i)=fi(1);
    if length(fi)>1, Mesh.boundright(i)=fi(2); end
end
Mesh.boundmarker(Mesh.boundright==0)=-1;

% tripatchmod(Mesh,rand(Mesh.ncells,1))
% if nargout<2, tripatchmod(Mesh); end