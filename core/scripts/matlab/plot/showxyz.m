function showxyz(P,mark)

% SHOWXYZ - show coordinate in matrix as points
% showxyz(P[,mark]);

if nargin<2, mark='.'; end
plot3(P(:,1),P(:,2),P(:,3),mark);
axis equal tight
grid on