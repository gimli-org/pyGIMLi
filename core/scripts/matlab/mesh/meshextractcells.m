function Mout=meshextractcells(Mesh,fi,delunused)

% MESHEXTRACTCELLS - Extract mesh cells
% Mesh = meshextractcells(Mesh,indices[,removeunusednodes])

if nargin<3, delunused=0; end
Mout=Mesh;
if islogical(fi), fi=find(fi); end
Mout.ncells=length(fi);
Mout.cell=Mesh.cell(fi,:);
fn=fieldnames(Mesh);
for i=1:length(fn),
  ff=getfield(Mesh,fn{i});
  if isequal(sort(size(ff)),[1 Mesh.ncells]), 
      Mout=setfield(Mout,fn{i},ff(fi)); end
end
Mout.ncells=size(Mout.cell,1);
if delunused,
    [un,uni]=unique(Mout.cell(:));
    [sun,loc]=ismember(un,1:Mesh.nnodes);
    Mout.nnodes=length(loc);
    map=zeros(Mesh.nnodes,1);
    map(loc)=1:length(loc);
    Mout.node=Mesh.node(loc,:);
    Mout.nodemarker=Mesh.nodemarker(loc);
    Mout.cell=map(Mout.cell);
end