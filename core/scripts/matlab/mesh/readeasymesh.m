function Mesh = readeasymesh(base)

% READEASYMESH - Reads name from easymesh output (.n/.e)
% Mesh = readeasymesh(basename)

if nargin<1, error('Specify mesh basename'); end
Mesh=[];Mesh.dim=2;
fid=fopen([base '.n'],'r');
if fid<0, error(['File ' base '.n does not exist']); end
Mesh.nnodes=fscanf(fid,'%d\n',1);
Mesh.node=zeros(Mesh.nnodes,Mesh.dim);
Mesh.nodemarker=zeros(Mesh.nnodes,1);
for i=1:Mesh.nnodes,
   zeile=strrep(fgetl(fid),':','');
   aa=str2num(zeile);
   Mesh.node(i,1:2)=aa(2:3);
   Mesh.nodemarker(i)=aa(4);
end
fclose(fid);
fid=fopen([base '.e'],'r');
if fid<0, error(['File ' base '.e does not exist']); end
Mesh.ncells=fscanf(fid,'%d\n',1);
Mesh.cell=zeros(Mesh.ncells,Mesh.dim+1);
Mesh.cellmarker=zeros(Mesh.ncells,1);
for i=1:Mesh.ncells,
   zeile=strrep(fgetl(fid),':','');
   aa=str2num(zeile);
   Mesh.cell(i,1:3)=aa(2:4);
   Mesh.cellmarker(i)=aa(12);
end
fclose(fid);
Mesh.cell=Mesh.cell+1;
fid=fopen([base '.s'],'r');
if fid<0, error(['File ' base '.s does not exist']); end
Mesh.nbounds=fscanf(fid,'%d\n',1);
Mesh.bounds=zeros(Mesh.nbounds,Mesh.dim);
Mesh.boundleft=zeros(Mesh.nbounds,1);
Mesh.boundright=Mesh.boundleft;
Mesh.boundmarker=Mesh.boundleft;
for i=1:Mesh.nbounds,
   zeile=strrep(fgetl(fid),':','');
   aa=str2num(zeile);
   Mesh.bound(i,1:2)=aa(2:3);
   Mesh.boundleft(i)=aa(4);
   Mesh.boundright(i)=aa(5);
   Mesh.boundmarker(i)=aa(6);   
end
fclose(fid);
Mesh.bound=Mesh.bound+1;
Mesh.boundleft=Mesh.boundleft+1;
Mesh.boundright=Mesh.boundright+1;
fprintf('%dD Mesh loaded (%d nodes, %d cells)\n',Mesh.dim,Mesh.nnodes,Mesh.ncells);

% clf;patch('Vertices',Mesh.node,'Faces',Mesh.cell,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]);
