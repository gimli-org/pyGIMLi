function Mesh = hex2tet(Mod,y,z)

% HEX2TET - Convert 3d grid model to 
% Mesh = hex2tet(Mod)
% Mod - 3d grid model containing x,y,z and M OR
% Mesh = hex2tet(x,y,z)

% Mesh - 3d tetraedral model containing node,cell and cellattr

Mesh=[];Mesh.dim=3;
if nargin>2,
    [X,Y,Z]=ndgrid(Mod,y,z);
    nx=length(Mod);ny=length(y);nz=length(z);
elseif isstruct(Mod),
    if ~(isfield(Mod,'x')&isfield(Mod,'y')&isfield(Mod,'z')),
        error('Model must contain x, y and z!');
    end
    [X,Y,Z]=ndgrid(Mod.x,Mod.y,Mod.z);
    nx=length(Mod.x);ny=length(Mod.y);nz=length(Mod.z);
end
Mesh.node=[X(:) Y(:) Z(:)];
Mesh.cell=[];
nxy=nx*ny;
% uninode=[1 2 7 3;1 3 7 4;1 7 8 4;6 7 1 2;6 5 1 7;5 8 1 7]; % 6 tetrahedra
uninode=[1 2 4 5;2 3 4 7;2 4 5 7;6 5 7 2;5 8 7 4]; % 5 tetrahedra
% uninode=uninode(:,[1 2 4 3]);
nt=size(uninode,1);
if isfield(Mod,'M'),
    aa=repmat(Mod.M(:)',nt,1);
    Mesh.cellattr=aa(:);
end
nm=(nx-1)*(ny-1)*(nz-1);
Mesh.cell=zeros(nm*size(uninode,1),4);
for k=1:nz-1,
    for j=1:ny-1,
        for i=1:nx-1,
            cellxy=[i i+1 i+nx+1 i+nx ]+ny*(j-1);
            cellp=[nxy*(k-1)+cellxy nxy*k+cellxy];
            nc=i-1+(nx-1)*(j-1)+(nx-1)*(ny-1)*(k-1);
            Mesh.cell(nc*nt+(1:nt),:)=cellp(uninode);
        end
    end
end
Mesh.nnodes=size(Mesh.node,1);
Mesh.ncells=size(Mesh.cell,1);
Mesh.nbounds=0;
Mesh.nodemarker=zeros(Mesh.nnodes,1);
Mesh.cellnodes=ones(Mesh.ncells,1)*4;