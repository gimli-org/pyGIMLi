function showpoly2d(Poly)

% SHOWPOLY2D - Show 2d poly file
% showpoly2d(Poly) shows nodes, edges and region marker
% edge colors: green (outside), red (surface), magenta (interior) & blue

cols={'green','red','blue','magenta'};
if isfield(Poly,'nodemarker'), fi=find(Poly.nodemarker==-99);
else fi=find(Poly.node(:,3)==-99); end
plot(Poly.node(:,1),Poly.node(:,2),'bx',Poly.node(fi,1),Poly.node(fi,2),'ro');
for i=1:size(Poly.edge,1),
    line(Poly.node(Poly.edge(i,1:2),1),Poly.node(Poly.edge(i,1:2),2),'Color',cols{min(Poly.edge(i,3)+3,length(cols))});
end
if isfield(Poly,'region'), 
    hold on;
    for i=1:size(Poly.region,1),
        plot(Poly.region(i,1),Poly.region(i,2),'x','MarkerSize',5);
        text(Poly.region(i,1),Poly.region(i,2),num2str(Poly.region(i,3)));
    end
    hold off;
end
axis equal tight xy