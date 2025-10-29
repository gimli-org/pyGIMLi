function vind = getcindex(val,minval,maxval,nind)

% GETCINDEX - Get color index for graphics patching
% vind = getcindex(values,minvalue,maxvalue[,ncolors])
% values .. vector of values
% minvalue/maxvalue .. the ends of the color scale
% ncolors .. number of colors for colormap (default 64)

if nargin<4, nind=64; end
vind=round((val-minval)/(maxval-minval)*nind);
vind=min(max(1,vind),nind);