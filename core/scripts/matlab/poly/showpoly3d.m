function showpoly3d(Poly)

% SHOWPOLY3d - Show 3d poly file
% showpoly3d(Poly)

cla;
for i=1:length(Poly.face),
    ff=Poly.face{i};
%     line(Poly.node(ff,1),Poly.node(ff,2),Poly.node(ff,3),'Color','red');
    patch(Poly.node(ff,1),Poly.node(ff,2),Poly.node(ff,3),[rand(1) 0.5 0.5],'EdgeColor',[1 0 0]);
end
hold on
showxyz(Poly.node,'k.')
if size(Poly.node,2)>=4, showxyz(Poly.node(Poly.node(:,4)==-99,:),'b.'); end
if isfield(Poly,'hole')&&~isempty(Poly.hole),
    showxyz(Poly.hole,'go');
end
if isfield(Poly,'region')~isempty(Poly.region),
    showxyz(Poly.region,'gx');
end
hold off
xlabel('x in m');ylabel('y in m');zlabel('z in m');