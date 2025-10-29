function savepointvtk(filename,points)

% SAVEPOINTVTK - Save vtk of points (e.g. electrodes)
%                for use in VTK rendering (e.g. paraview)
% savepointvtk(filename,points)

if nargin<2, error('Specify filename and point array'); end
if size(points,2)<3, points(1,3)=0; end

fid=fopen(filename,'w');
fprintf(fid,'# vtk DataFile Version 3.0\n');
fprintf(fid,'3d point set\nASCII\nDATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'POINTS %d double\n',size(points,1));
fprintf(fid,'%f\t%f\t%f\n',points');
fclose(fid);