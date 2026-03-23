function writeconfigfile(filename,AA,header)

% WRITECONFIGFILE = write config file from structure variable
% writeconfigfile(filename,STRUC[,headerstring])
% Example: STRUC=struct('a',1,'b',23.1)
% STRUC =
%   a: 1
%   b: 23.1
% Then the file looks like
% a=1
% b=23.1
% See also readconfigfile

if ~isstruct(AA), return; end
fid=fopen(filename,'w');
if fid==-1, display(['error opening file ' filename]);return; end
if (nargin>2)&&ischar(header), fprintf(fid,'#%s\n',header); end
fn=fieldnames(AA);
for i=1:length(fn),
    ff=getfield(AA,fn{i});
    if isnumeric(ff), fprintf(fid,'%s=%g\n',fn{i},ff); 
    elseif ischar(ff), fprintf(fid,'%s=%s\n',fn{i},ff); 
    end
end
fclose(fid);
