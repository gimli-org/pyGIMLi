function Mesh = grid4mesh(x,z,M)

% GRID2MESH - Transform regular (2d) grid to triangle mesh
%             by dividing each rectangle into 4 triangles
% [Mesh,T] = grid2mesh(x,z[,M])

if nargin<2, error('Too less input arguments (x,y)!'); end
if nargin<3, M=ones(length(x)-1,length(z)-1); end

Mesh=[];Mesh.dim=2;
nx=length(x);nz=length(z);
x=x(:);z=z(:);
dx=diff(x);dz=diff(z);
Mesh.nnodes=nx*nz;
Mesh.node=zeros(Mesh.nnodes,Mesh.dim);
for i=1:nz,
    Mesh.node((i-1)*nx+1:nx*i,1)=x(:);
    Mesh.node((i-1)*nx+1:nx*i,2)=z(i);
end
Mesh.nodemarker=zeros(Mesh.nnodes,1);
Mesh.ncells=(nx-1)*(nz-1);
Mesh.cell=zeros(Mesh.ncells,Mesh.dim*2);
Mesh.cellnodes=ones(Mesh.ncells,1)*(Mesh.dim*2);
Mesh.nbounds=(nx-1)*nz+(nz-1)*nx;
Mesh.boundnodes=ones(Mesh.nbounds,1)*2;
Mesh.boundmarker=zeros(Mesh.nbounds,1);
Mesh.bound=zeros(Mesh.nbounds,Mesh.dim);
Mesh.bound(1:nx-1,1)=(1:nx-1)';Mesh.bound(1:nx-1,2)=(2:nx)';
ic=0;ib=1;
iic=[0 1 nx+1 nx];
iib=[0 1;0 nx];
for iz=1:nz-1,
    for ix=1:nx-1,
        ic=ic+1;
        offset=ix+(iz-1)*nx;
        Mesh.cell(ic,:)=iic+offset;
        Mesh.bound(ib:ib+1,1:2)=iib+offset;
        ib=ib+2;
    end
    Mesh.bound(ib,1:2)=[1 nx+1]+offset;
    ib=ib+1;
end
li=(1:nx-1)'+nx*(nz-1);
Mesh.bound(ib:ib+nx-2,1:2)=[li li+1];

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
if nargin>3, Mesh.cellattr=M(:); else Mesh.cellattr=ones(Mesh.ncells,1)*2; end

% tripatchmod(Mesh,rand(Mesh.ncells,1))
if nargout<1, tripatchmod(Mesh); end