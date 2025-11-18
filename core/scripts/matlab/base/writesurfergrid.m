function loadsurfergrid(fname,A,x,y);

% LOADSURFERGRID - Load surfer grid (*.grd) file
% [A,x,y] = loadsurfergrid(filename);

if nargin<2, error('specify matrix!'); end
if nargin<4, y=size(A,2); end
if nargin<3, x=size(A,2); end

fid=fopen(fname,'w');
if fid<0, display('Could not open file'); return; end
nbytes1=4;
nbytes2=72;
vers=1;
xll=x(1);
ncol=length(x);
xsize=rndig(median(diff(x)));
yll=y(1);
nrow=length(y);
ysize=rndig(median(diff(y)));
zmin=min(A(:));
zmax=max(A(:));
rotation=0;
blankvalue=0.0;
eins=77568;
% Header section
fwrite(fid,'DSRB','char');
fwrite(fid,nbytes1,'long');
fwrite(fid,vers,'long');
% Grid section
fwrite(fid,'GRID','char');
fwrite(fid,nbytes2,'long');
fwrite(fid,nrow,'long');
fwrite(fid,ncol,'long');
fwrite(fid,xll,'double');
fwrite(fid,yll,'double');
fwrite(fid,xsize,'double');
fwrite(fid,ysize,'double');
fwrite(fid,zmin,'double');
fwrite(fid,zmax,'double');
fwrite(fid,rotation,'double');
fwrite(fid,blankvalue,'double');
% Data section
fwrite(fid,'DATA','char');
fwrite(fid,eins,'long')
fwrite(fid,A','double')';
fclose(fid);