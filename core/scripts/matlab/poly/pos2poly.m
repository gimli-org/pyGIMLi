function Poly=pos2poly(pos,maxdep,dd)

% POS2POLY - Creates poly file from position matrix
% Poly = pos2poly(pos,maxdep,dx)

if nargin<3, dd=0.2; end
if nargin<2, maxdep=(pos(end,1)-pos(1,1))/4; end

xx=pos(:,1)';
if dd>0,
    A=[xx-[0 dd*diff(xx)];xx;xx+[dd*diff(xx) 0]];
    if dd>0.333, A(1,:)=[]; end
    xx=A(:);
else
    xx=[xx(1);xx(:);xx(end)];
end
zz=interp1(pos(:,1),pos(:,2),xx);
if dd>0.333, xx=[xx(1);xx];zz=[zz(1);zz]; end
zz(1)=zz(1)-maxdep;zz(end)=zz(end)-maxdep;
Poly.node=[xx zz];
if dd==0,
    Poly.node(2:end-1,3)=-99;
elseif dd>0.333, 
    Poly.node(2:2:end-1,3)=-99;
else,
    Poly.node(2:3:end-1,3)=-99; 
end
% Poly.node([1 end],3)=0;
Poly.edge=(1:size(Poly.node,1))';
Poly.edge(:,2)=Poly.edge(:,1)+1;Poly.edge(end,2)=1;
Poly.edge(:,3)=0;