function writepoly2d(filename,Poly)

% WRITEPOLY2d - Create poly file from Poly struct
% writepoly2d(filename,Poly)
% Poly.node .. nodes with x y and node marker
%     (-99 electrode, -999/-1000 reference current/potential, 0 else)
% Poly.edge .. edges with node1, node2 and edge marker
%     (-1 Neumann, -2 Dirichlet, >1 user line, 0 else = Neumann)
% Poly.region .. region marker with x, y, marker and maxsize)

np=size(Poly.node,1);
nn=zeros(np,4);
nn(:,1)=(1:np)'-1;
nn(:,2:size(Poly.node,2)+1)=Poly.node;
ne=size(Poly.edge,1);
ee=zeros(ne,4);
ee(:,1)=(1:ne)'-1;
ee(:,2:size(Poly.edge,2)+1)=Poly.edge-1;
ee(:,end)=ee(:,end)+1;
fid=fopen(filename,'w');
fprintf(fid,'%d %d %d %d\n',np,2,0,1);
fprintf(fid,'%d\t%.3f\t%.3f\t%d\n',nn');
fprintf(fid,'%d %d\n',ne,1);
fprintf(fid,'%d\t%d\t%d\t%d\n',ee');
if isfield(Poly,'region'),
    fprintf(fid,'0\n%d\n',size(Poly.region,1));
    if size(Poly.region,2)>3,
        for i=1:size(Poly.region),
            fprintf(fid,'%d %g %g %g %g\n',i-1,Poly.region(i,:));
        end        
    else
        for i=1:size(Poly.region),
            fprintf(fid,'%d %g %g %g %g\n',i-1,Poly.region(i,:),0);
        end
    end
else    
    fprintf(fid,'%d\n%d\n',0,0);
end
fclose(fid);
% clf;
% plot(nn(:,2),nn(:,3),'*');
% for i=1:size(ee,1),
%     line(nn(ee(i,2:3)+1,2),nn(ee(i,2:3)+1,3));
% end
