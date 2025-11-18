function [Mesh,Poly]=createparamesh2d(pos,basename,quality,pdep,pbou,bou,dep,dd)

% CREATEPARAMESH - Create para mesh (outer space=1,para=2) from positions
% createparamesh2d(pos,basename,pdep,pbou,bou,dep);
% pdep..paradepth
% pbou..paraboundary (in m (>0), in % (<0), automatic (==0)
% bou..boundary on the whole domain
% dep..depth of the whole domain

if nargin<2, basename='mesh'; end
if nargin<3, quality=33.8; end
if nargin<4, pdep=0; end
if nargin<5, pbou=0; end
if nargin<6, bou=0; end
if nargin<7, dep=bou; end
if nargin<8, dd=1; end

%automatic
ext=(max(pos(:,1))-min(pos(:,1)));
if pdep==0, pdep=ext/4; end
if pbou==0, pbou=ext*0.1; end
if pbou<0, pbou=exp*abs(pbou)/100; end

%% create nodes
Poly.node=pos;Poly.node(:,3)=-99;
di=sqrt(sum(diff(pos).^2,2));
if dd>0, 
    fi=find(di<dd*2);
    newnode=pos(1:end-1,1)+dd;newnode(:,3)=0;
    newnode(:,2)=interp1(pos(:,1),pos(:,2),newnode(:,1));
    Poly.node=sortrows([Poly.node;newnode]);
end
Poly.node(end+1,1:2)=Poly.node(end,1:2)+[pbou 0];
Poly.node(end,2)=interp1(pos(:,1),pos(:,2),Poly.node(end,1),'linear','extrap');
Poly.node(end+1,1:2)=Poly.node(end,1:2)-[0 pdep];
Poly.node(end+1,1:2)=Poly.node(1,1:2)-[pbou pdep];
Poly.node(end+1,1:2)=Poly.node(1,1:2)-[pbou 0];
Poly.node(end,2)=interp1(pos(:,1),pos(:,2),Poly.node(end,1),'linear','extrap');
Poly.edge=(1:size(Poly.node,1))';
Poly.edge(:,2)=Poly.edge(:,1)+1;Poly.edge(end,2)=1;
nn=size(Poly.node,1);
Poly.node(end+1,1:2)=Poly.node(nn-3,1:2)+[bou 0];
Poly.node(end,2)=interp1(pos(:,1),pos(:,2),Poly.node(end,1),'linear','extrap');
Poly.node(end+1,1:2)=Poly.node(end,1:2)-[0 dep];
Poly.node(end+1,1:2)=Poly.node(nn,1:2)-[bou dep];
Poly.node(end+1,1:2)=Poly.node(nn,1:2)-[bou 0];
Poly.node(end,2)=interp1(pos(:,1),pos(:,2),Poly.node(end,1),'linear','extrap');
Poly.node(end-1,2)=Poly.node(end,2)-dep;
%% create edged
Poly.edge=[Poly.edge;nn-3 nn+1;nn+1 nn+2;nn+2 nn+3;nn+3 nn+4;nn+4 nn];
Poly.region=[Poly.node(end,1:2)+[1 -1] 1;Poly.node(1,1:2)-[0 0.1] 2];
writepoly2d([basename '.poly'],Poly)
dos(['dctriangle -v -q' num2str(quality) ' ' basename '.poly']);
Mesh=loadmesh([basename '.bms']);
fprintf('Found mesh with %d parameters (total %d cells)',...
    length(find(Mesh.cellattr==2)),Mesh.ncells);