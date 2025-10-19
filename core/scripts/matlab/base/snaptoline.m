function [l,x1,y1] = snapline(x,y)

% SNAPLINE - SNAPS points (x,y) onto a straight line
% [l,x1,y1] = snapline(x,y)

if size(x)~=size(y), error('size does not match!'); end
nl=length(x);
G=ones(nl,2);
G(:,1)=x;
ab=G\y,a=ab(1);b=ab(2);
%% description by parametric form
%% (x,y) = (0,b) + (1,a) * t
x1=(x+y*a-a*b)/(a*a+1);
y1=x1*a+b;
xx=[min(x) max(x)];yy=xx*a+b;
plot(x,y,'bx',x1,y1,'ro',xx,yy,'r-');axis equal
dr=sqrt((x1-x).^2+(y1-y).^2);
%% profile length
l=(x-x(1))*(1+a^2);
if max(l)<=0, l=-l; end
