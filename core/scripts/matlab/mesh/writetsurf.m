function writetsurf(filename,Mesh,xy,varname)

% WRITETSURF - Write TSurf Mesh to file (*.ts)
% writetsurf(filename,Mesh)
%   Mesh.patt contains attributes for each point/node
%   alternatively Mesh.att (attributes for each cell) is used
% writetsurf(filename,Mesh,xy) with xy=point list for 2d->3d
% writetsurf(filename,Mesh,xy,varname) with variable name string

if nargin<1, filename='test.ts'; end
if nargin<2, error('Specify mesh!'); end
if nargin<3, xy=[0 0;1 0]; end
if nargin<4, varname='resistivity'; end

if ~isfield(Mesh,'patt')&&isfield(Mesh,'cellattr'), % interpolate
    Mesh.patt=zeros(Mesh.nnodes,1);
    for i=1:Mesh.nnodes,
       [row,col]=find(Mesh.cell==i);
       Mesh.patt(i)=mean(Mesh.cellattr(row));
    end
end
if ~isfield(Mesh,'patt')&&isfield(Mesh,'att'), % interpolate
    Mesh.patt=zeros(Mesh.nnodes,1);
    for i=1:Mesh.nnodes,
       [row,col]=find(Mesh.cell==i);
       Mesh.patt(i)=mean(Mesh.att(row));
    end
end


name=strrep(filename,'.ts',''); % comment

node3d=Mesh.node;
if size(Mesh.node,2)==2, 
    node3d(:,3)=Mesh.node(:,2);
    node3d(:,2)=0;
end
if nargin>2, % xy points given
    mbm=[0;cumsum(sqrt(sum(diff(xy).^2,2)))];
    node3d(:,1)=interp1(mbm,xy(:,1),Mesh.node(:,1),'linear','extrap');
    node3d(:,2)=interp1(mbm,xy(:,2),Mesh.node(:,1),'linear','extrap');
end

fid=fopen(filename,'w');
fprintf(fid,'GOCAD TSURF 1\r\nHEADER {\r\n');
fprintf(fid,'name:%s\r\n',name);
fprintf(fid,'mesh:on\r\n*solid*color:1 0.447059 0.337255 1\r\n');
fprintf(fid,'ivolmap:false\r\nimap:false\r\n*painted:on\r\n');
fprintf(fid,'*painted*variable:%s\r\nlast_selected_folder:Texture\r\n',varname);
fprintf(fid,'}\r\n');
fprintf(fid,'PROPERTIES %s\r\nPROP_LEGAL_RANGES **none**  **none**\r\n',varname);
fprintf(fid,'NO_DATA_VALUES -99999\r\nPROPERTY_CLASSES test\r\nPROPERTY_KINDS Acoustic Impedance\r\n');
fprintf(fid,'PROPERTY_SUBCLASSES QUANTITY Float\r\nESIZES 1\r\nUNITS Ohmm');
fprintf(fid,'PROPERTY_CLASS_HEADER test {\r\n*low_clip:%f\r\n',interperc(Mesh.patt,5));
fprintf(fid,'*high_clip:%f\r\n*pclip:99\r\n}\r\n',interperc(Mesh.patt,95));

fprintf(fid,'TFACE\r\n');
fmt='PVRTX %d %f %f %f';
mat=[(1:Mesh.nnodes);node3d'];
if isfield(Mesh,'patt'), fmt=[fmt ' %f'];mat=[mat;Mesh.patt(:)']; end
fprintf(fid,[fmt '\r\n'],mat);
fprintf(fid,'TRGL %d %d %d\r\n',Mesh.cell');
% fprintf(fid,'BSTONE 1\r\n');
fprintf(fid,'END');
fclose(fid);