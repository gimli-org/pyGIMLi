function c = b2r(m)
% Create Blue-white-red colormap(for zero-symmetric data)
% c=b2r(m)   m..size of colormap (default from current figure)

if nargin < 1, m = size(get(gcf,'colormap'),1); end
n = fix(m/2);

r = [sqrt(1:n-1)'/sqrt(n-1); ones(m-n+1,1)];
g = [sqrt(1:n-1)'/sqrt(n-1);  ones(1,1); flipud(sqrt(1:n)'/sqrt(n))];
b = flipud(r);
length(r);
length(g);
length(b);


c = [r g b];
