function cellmids=meshcellmid(Mesh)

% MESHCELLMID - Return cell midpoint of mesh
% cellmids = meshcellmid(Mesh)
% Mesh consists of nodes (n x dim) and cells (c x nodepercell)

cellmids=zeros(Mesh.ncells,size(Mesh.node,2));
for i=1:Mesh.ncells,
    cellmids(i,:)=mean(Mesh.node(Mesh.cell(i,:),:));
end