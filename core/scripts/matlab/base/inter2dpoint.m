function poz = inter2dpoint(to,po)

% INTER2DPOINT - Interpolate 2d to points by delaunay linear interpolation
% z = inter2dpoint(xyz,xy)

tri=delaunay(to(:,1),to(:,2));
t=tsearch(to(:,1),to(:,2),tri,po(:,1),po(:,2));
poz=zeros(size(po,1),1);
for i=1:size(po,1),
    x=to(tri(t(i),:),1);y=to(tri(t(i),:),2);z=to(tri(t(i),:),3);
    x21=x(2)-x(1);y21=y(2)-y(1);x31=x(3)-x(1);y31=y(3)-y(1);
    jdet=x21*y31-x31*y21;
    xp1=po(i,1)-x(1);yp1=po(i,2)-y(1);
    r = ( xp1 * y31 - x31 * yp1 ) / jdet;
    s = ( x21 * yp1 - xp1 * y21 ) / jdet;
    poz(i)=(1-r-s)*z(1)+r*z(2)+s*z(3);    
end