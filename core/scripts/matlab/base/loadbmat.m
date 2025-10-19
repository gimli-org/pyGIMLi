function S = loadbmat(sensname,inttype)

% LOADSENS - Load binary matrix file
% S = loadbmat(filename)

if nargin<1, error('specify filename'); end
if nargin<2, inttype='long'; end
fid=fopen(sensname,'r');
ndata=fread(fid,1,inttype);
nmodel=fread(fid,1,inttype);
S=zeros(ndata,nmodel);
for i=1:ndata, S(i,:)=fread(fid,nmodel,'double')'; end
fclose(fid);