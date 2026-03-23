function savematrixbinary(filename,A)

% SAVEMATRIXBINARY - Save matrix to binary (bmat) format
% savesensmat(S,sensname)

if nargin<2, error('Two inputs (name and matrix) required!'); end
fid=fopen(filename,'w');
ndata=size(A,1);nmodel=size(A,2);
fwrite(fid,ndata,'long');
fwrite(fid,nmodel,'long');
for i=1:ndata, fwrite(fid,A(i,:),'double'); end
fclose(fid);
