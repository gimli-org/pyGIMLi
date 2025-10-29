function Mesh=loadmesh(meshname,dobound)

% LOADMESH - Loads DCFEMLIB mesh from file(s)
%            mesh may be 2d/3d and binary or ascii
% Mesh = loadmesh(meshname);
% Mesh = loadmesh(meshname[,dobound]); %calculates boundaries new

if nargin<2, dobound=0; end
Mesh=[];if nargin<1, meshname='testmesh2_h2.bms'; end
if exist([meshname '.bms'],'file'),
    meshname=[meshname '.bms'];
elseif exist([meshname '.n'],'file'),
%     fid=fopen([meshname '.n']);aa=str2num(fgetl(fid));fclose(fid);
%     if length(aa)==2, 
    [Mesh.cell,Mesh.node]=loadmesh2d(meshname);
    if size(Mesh.cell,2)>3, Mesh.cell(:,4:end)=[]; end
    Mesh.dim=2;
    Mesh.ncells=size(Mesh.cell,1);
    Mesh.nnodes=size(Mesh.node,1);
    Mesh.nbounds=0;
    Mesh.cellnodes=ones(Mesh.ncells,1)*3;
    Mesh.cellmarker=zeros(Mesh.ncells,1);
    Mesh.nodemarker=zeros(Mesh.nnodes,1);
    Mesh.cellattr=(1:Mesh.ncells)';
    return;
end

if ~exist(meshname,'file'),
    Mesh=[];
    display(['Filename ' meshname ' not found!']);
    return;
end
zahl='int32';wert='double';
fid=fopen(meshname);
if fid<0, error('Could not open file!'); end
Mesh.dim=fread(fid,1,zahl);
vertinfo=fread(fid,127,zahl);
Mesh.nnodes=fread(fid,1,zahl);
Mesh.node=fread(fid,[Mesh.dim Mesh.nnodes],wert)';
Mesh.nodemarker=fread(fid,Mesh.nnodes,zahl);
cellinfo=fread(fid,127,zahl);
Mesh.ncells=fread(fid,1,zahl);
Mesh.cellnodes=fread(fid,Mesh.ncells,zahl);
Mesh.cell=ones(Mesh.ncells,max(Mesh.cellnodes));
for i=1:Mesh.ncells, 
    nn=Mesh.cellnodes(i);
    Mesh.cell(i,1:nn)=fread(fid,nn,zahl)+1; % 1 counter
end
Mesh.cellattr=fread(fid,Mesh.ncells,wert);
boundinfo=fread(fid,127,zahl);
Mesh.nbounds=fread(fid,1,zahl);
if Mesh.nbounds>0,
    Mesh.boundnodes=fread(fid,Mesh.nbounds,zahl);
    Mesh.bound=ones(Mesh.nbounds,max(Mesh.boundnodes));
    for i=1:Mesh.nbounds,
        nn=Mesh.boundnodes(i);
        Mesh.bound(i,1:nn)=fread(fid,nn,zahl)+1;
    end
    Mesh.boundmarker=fread(fid,Mesh.nbounds,zahl);
    Mesh.boundleft=fread(fid,Mesh.nbounds,zahl)+1;
    Mesh.boundright=fread(fid,Mesh.nbounds,zahl)+1;
end
fclose(fid);
fprintf('%dD Mesh loaded (%d nodes, %d cells)\n',Mesh.dim,Mesh.nnodes,Mesh.ncells);
if dobound==0, return; end
if ~isfield(Mesh,'bound'),
    if Mesh.dim==3,
        Mesh.nbounds=Mesh.ncells*4;
        Mesh.bound=zeros(Mesh.nbounds,3);
        for i=1:Mesh.ncells,
            Mesh.bound((i-1)*4+(1:4),:)=[Mesh.cell(i,[1 1 1 2])'...
                                         Mesh.cell(i,[2 2 3 3])'...
                                         Mesh.cell(i,[3 4 4 4])'];
        end
    end
    if Mesh.dim==2,
        Mesh.nbounds=Mesh.ncells*3;
        Mesh.bound=zeros(Mesh.nbounds,2);
        for i=1:Mesh.ncells,
            Mesh.bound((i-1)*3+(1:3),:)=[Mesh.cell(i,[1 1 2])' Mesh.cell(i,[2 3 3])'];
        end
    end
end
if (Mesh.dim==3), % edges for 3d mesh (Dijkstra algorithm
    Mesh.edge=unique(sort([Mesh.bound(:,1:2);Mesh.bound(:,1:2:3);Mesh.bound(:,2:3)],2),'rows');
    Mesh.edgelength=sqrt(sum((Mesh.node(Mesh.edge(:,1),:)-Mesh.node(Mesh.edge(:,2),:)).^2,2));
    Mesh.edgeslowness=ones(size(Mesh.edgelength));
    Mesh.edgeneigh=zeros(length(Mesh.edge),5);
    for i=1:length(Mesh.edge), % hier ist noch einiges kaputt
        [i1,j1]=find(Mesh.cell==Mesh.edge(i,1));
        [i2,j2]=find(Mesh.cell==Mesh.edge(i,2));
        fi=intersect(i1,i2);
        Mesh.edgeneigh(i,1:length(fi))=fi;
    end

%     Mesh.boundneigh=zeros(Mesh.nbounds,4);
%     for i=1:Mesh.nbounds,
%         [i1,j1]=find(Mesh.cell==Mesh.bound(i,1));
%         [i2,j2]=find(Mesh.cell==Mesh.bound(i,2));
%         ii=intersect(i1,i2);
%         Mesh.boundneigh(i,1:length(ii))=ii;
%     end
%     [mm,ii]=unique(Mesh.bound,'rows');
%     Mesh.bound=mm;Mesh.boundneigh=Mesh.boundneigh(ii,:);
%     Mesh.nbounds=length(ii);
end
if(~isfield(Mesh,'boundleft'))||(max(Mesh.boundleft)<=0), % find left/right elements to edges
    Mesh.boundleft=zeros(Mesh.nbounds,1);
    Mesh.boundright=zeros(Mesh.nbounds,1);
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
end
return
cla reset;
for i=1:Mesh.ncells,
    nums=Mesh.cell(i,:);
    patch(Mesh.node(nums,1),Mesh.node(nums,2),[1 1 1]*(i-1)/Mesh.ncells);
end
% //   long[ 1 ] dimension
% //   long[ 127 ] dummy vertices information
% //   long[ 1 ] nVerts, number of vertices
% //   double[ dimension * nVerts ]; koordinates,  dimension == 2 ( x, y ), dimension == 3 (x, y, z )
% //   long[ nVerts ] vertex markers
% //   long[ 127 ] dummy cell information
% //   long[ 1 ] nCells, number of cell
% //   long[ nCells ] cellVerts; nodes for each cell
% //   long[ sum( cellVerts ) ] cellsidx
% //   double[ nCells ] attribute, cell attributes
% //   long[ 127 ] dummy boundary information
% //   long[ 1 ] nBounds, number of boundarys
% //   long[ nBounds ] boundVerts ; nodes for each boundary
% //   long[ sum( boundVerts ) ] boundIdx
% //   long[ nBounds ] boundary markers
% //   long[ nBounds ] leftNeighbour idx (-1) if no neighbour present or info unavailable
% //   long[ nBounds ] rightNeighbour idx (-1) if no neighbour present or info unavailable