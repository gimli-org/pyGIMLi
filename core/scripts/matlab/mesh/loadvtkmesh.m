function Mesh=loadvtkmesh(filename,dim)

% LOADVTKMESH - Load mesh from VTK format
% loadvtkmesh(filename)
% loadvtkmesh(filename,dim)

if nargin<2, dim=2; end
Mesh=[];
if nargin<1, filename='showmodel.vtk'; end
fid=fopen(filename,'r');
for i=1:5, zeile=fgetl(fid); end
Mesh.nnodes=sscanf(zeile,'%*s%d*%s');
Mesh.node=fscanf(fid,'%f',[3 Mesh.nnodes])';
zeile='';
while length(zeile)<4, zeile=fgetl(fid); end
Mesh.ncells=sscanf(zeile,'%*s%d*%s');
if nargin>1, % dim given
    frow=dim+2;
else
    fpos=ftell(fid);
    zeile=fgetl(fid);
    nn=str2num(zeile);
    frow=length(nn);
    dim=frow-2;
    fseek(fid,fpos,'bof');
end
Mesh.cell=fscanf(fid,'%f',[frow Mesh.ncells])'+1;
Mesh.cell(:,1)=[]; % node number
zeile='';
while length(zeile)<2, zeile=fgetl(fid); end
nctypes=sscanf(zeile,'%*s%d');
ntypes=fscanf(fid,'%d',nctypes);
zeile='';
while length(zeile)<2, zeile=fgetl(fid); end
zeile='';while isempty(strfind(zeile,'LOOKUP')), zeile=fgetl(fid); end
Mesh.cellattr=fscanf(fid,'%f',Mesh.ncells);
zeile='';
while isstr(zeile)&&(length(zeile)<2), zeile=fgetl(fid); end
zeile=fgetl(fid);
Mesh.cellattr2=fscanf(fid,'%f',Mesh.ncells);
zeile='';
while isstr(zeile)&&(length(zeile)<2), zeile=fgetl(fid); end
zeile=fgetl(fid);
Mesh.cellattr3=fscanf(fid,'%f',Mesh.ncells);
zeile='';
while isstr(zeile)&&(length(zeile)<2), zeile=fgetl(fid); end
zeile=fgetl(fid);
Mesh.cellattr4=fscanf(fid,'%f',Mesh.ncells);
zeile='';
while isstr(zeile)&&(length(zeile)<2), zeile=fgetl(fid); end
zeile=fgetl(fid);
Mesh.cellattr5=fscanf(fid,'%f',Mesh.ncells);
fclose(fid);
Mesh.dim=dim;
Mesh.cellnodes=ones(Mesh.ncells,1)*(dim+1);