function saveunifile(fname,N)

% SAVEUNIFILE - Save unified data format
% saveunifile(filename,Data)

if nargin<2, error('saveunifile(filename,Data)'); end
if (nargin>1)&&isstruct(fname)&&ischar(N), %reversed order
    dummy=N;N=fname;fname=dummy;    
end
newline='\r\n';
fid=fopen(fname,'w');
if fid<0,
    error('File not found!');
end
%% positions
ss='#';
if isfield(N,'x')|isfield(N,'y'),
    elec=[];
    if isfield(N,'x'), elec=[elec N.x];ss=[ss '\tx']; end
    if isfield(N,'y'), elec=[elec N.y];ss=[ss '\ty']; end
    if isfield(N,'z'), elec=[elec N.z];ss=[ss '\tz']; end    
elseif isfield(N,'pos'),
    elec=N.pos;
else
    elec=N.elec;
end
fprintf(fid,'%d',size(elec,1));
fprintf(fid,'# Number of sensors');
fprintf(fid,newline);
fs='';
for i=1:size(elec,2), fs=[fs '\t%g']; end
ndata=0;
if length(ss)>1,
    ss(2:3)='';ss=[ss newline];
    fprintf(fid,ss);
end
fs(1:2)='';fs=[fs newline];
fprintf(fid,fs,elec');
%% data 
if isfield(N,'a'), ndata=length(N.a); end
if isfield(N,'s'), ndata=length(N.s); end
fprintf(fid,'%d',ndata);
fprintf(fid,'# Number of data');
fprintf(fid,newline);
ss='#';fs='';
fn=fieldnames(N);
DATA=[];
for i=1:length(fn),    
    ff=getfield(N,fn{i});
    if (min(size(ff))==1)&(length(ff)==ndata),
        DATA=[DATA ff(:)];
        fni=fn{i};
        if strcmp(fni,'r'), fni='rhoa'; end
        if strcmp(fni,'rho'), fni='R'; end
        if strcmp(fni,'imp'), fni='R'; end
        if strcmp(fni,'konf'), fni='k'; end
        ss=[ss '\t' fni];
        fs=[fs '\t%g'];
    end
end
ss(2:3)='';ss=[ss newline];
fs(1:2)='';fs=[fs newline];
fprintf(fid,ss);
fprintf(fid,fs,DATA');
if isfield(N,'topo'),
    fprintf(fid,['%d # topo points' newline],size(N.topo,1));
    tdim=size(N.topo,2);
    if tdim==2, 
        fprintf(fid,['#x z' newline]); 
        fprintf(fid,['%g\t%g' newline],N.topo');
    end
    if tdim==3, 
        fprintf(fid,['#x y z' newline]); 
        fprintf(fid,['%g\t%g\t%g' newline],N.topo');
    end    
end
fclose(fid);