function Poly=createcylpoly(rad,xli,nseg)

% CREATECYLPOLY - Create Cylinder domain poly file
% Poly = createcylpoly(rad,xlimits,nsegments)
%        both arguments can be neglected
%        rad..radius of cylinder
%        xlim..z-coordinates of top/bottom
%        nsegments..number of circle segments

if nargin<1, rad=1; end
if nargin<2, xli=[0 1]; end
if nargin<3, nseg=48; end

% rad=nseg=24;xli=[-0.8 0];

phi=(0:nseg-1)'/nseg*2*pi;
xy=[cos(phi) sin(phi)]*rad;
xy(:,3)=xli(1);
Poly.node=xy;
xy(:,3)=xli(2);
Poly.node=[Poly.node;xy];
for i=1:nseg,
    j=i+1;if j>nseg, j=1; end
    Poly.face{i}=[i j j+nseg i+nseg];
end
Poly.node(:,4)=0;
Poly.face{end+1}=1:nseg;
Poly.face{end+1}=(1:nseg)+nseg;
dz=abs(diff(xli))/4;
Poly.node(end+1,:)=[0 0 mean(xli)-dz/2 -999]; % reference electrode (arbitrary)
Poly.node(end+1,:)=[0 0 mean(xli)+dz/2 -1000]; % reference point for Neumann problem
if nargout<1, showpoly3d(Poly); end