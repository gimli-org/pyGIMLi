function meshareas=meshcellsize(Mesh)

% MESHCELLSIZE - Return mesh triangle areas or tetrahedra volume (dim=3)
% areas = meshcellsize(Mesh)
% Mesh consists of nodes (n x dim) and cells (c x nodepercell)

meshareas=zeros(Mesh.ncells,1);
if Mesh.dim==2,
    for i=1:Mesh.ncells,
        nodes=Mesh.node(Mesh.cell(i,:),:);
        a=sqrt(sum(diff(nodes(1:2,:)).^2));
        b=sqrt(sum(diff(nodes(2:3,:)).^2));
        c=sqrt(sum(diff(nodes(1:2:3,:)).^2));
        s=(a+b+c)/2;
        meshareas(i)=sqrt(s*(s-a)*(s-b)*(s-c));
    end
end
if Mesh.dim==3,
    for i=1:Mesh.ncells,
        nodes=Mesh.node(Mesh.cell(i,:),:);
        u2=sum((nodes(1,:)-nodes(2,:)).^2);
        U2=sum((nodes(3,:)-nodes(4,:)).^2);
        v2=sum((nodes(1,:)-nodes(3,:)).^2);
        V2=sum((nodes(2,:)-nodes(4,:)).^2);
        w2=sum((nodes(1,:)-nodes(4,:)).^2);
        W2=sum((nodes(2,:)-nodes(3,:)).^2);
        v=sqrt(4*u2*v2*w2 - u2*(v2+w2-U2)^2 - v2*(w2+u2-V2)^2 - w2*(u2+v2-W2)^2 +...
            (v2+w2-U2)*(w2+u2-V2)*(u2+v2-W2))/12;
        meshareas(i)=v;
    end
end