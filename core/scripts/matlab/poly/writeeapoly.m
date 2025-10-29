function writeeapoly(Poly,filename)

% WRITE Easymesh poly (*.d) file
% writeeapoly(Poly,filename)

if nargin<1, error('Specify poly file'); end
if nargin<2, filename='out.d'; end

%% write poly file
zerobased=1;
nnodes=size(Poly.node,1);
nedges=size(Poly.edge,1);
nn=ones(nnodes,5);
nn(:,1)=(1:nnodes)'-zerobased;
nn(:,2:size(Poly.node,2)+1)=Poly.node;
ee=ones(nedges,4);
ee(:,1)=(1:nedges)'-zerobased;
ee(:,2:size(Poly.edge,2)+1)=Poly.edge-zerobased;
[fid,mess]=fopen(filename,'w');
if fid>0,
    fprintf(fid,'%d\n',nnodes);
    fprintf(fid,'%d:\t%g\t%g\t%g\t%d\n',nn');
    fprintf(fid,'%d\n',nedges);
    fprintf(fid,'%d:\t%g\t%g\t%d\n',ee');
    fclose(fid);
else
    display(mess);
end