function R = turnmatrix(phixy,phiyz,phixz,isdegree)

% TURNMATRIX - Return matrix for koordinate system
% R = turnmatrix(phixy,phiyz,phixz[,isdegree])

if nargin<4, isdegree=-1; end
if nargin<3, phiyz=0; end
if nargin<2, phiyz=0; end
if isdegree==-1, isdegree=(max(abs([phixy phiyz phixz]))>pi); end
if isdegree,
   phixy=phixy*pi/180;
   phixz=phixz*pi/180;
   phiyz=phiyz*pi/180;
end

A=[cos(phixy) -sin(phixy) 0;sin(phixy) cos(phixy) 0;0 0 1];
B=[1 0 0;0 cos(phiyz) -sin(phiyz);0 sin(phiyz) cos(phiyz)];
C=[cos(phixz) 0 -sin(phixz);0 1 0;sin(phixz) 0 cos(phixz)];
R=A*B*C;