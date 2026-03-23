function dm=tantransdiff(p,l,u)

% TANTRANSDIFF - Derivative of tangens transformation with lower&upper bound
% m = tantransdiff(p,lowerbound,upperbound)
% both lower and upper bound are required
% See also TANTRANS, TANTRANSINV

if nargin<3, 
    error('Specify model, lower, and upper bound!'); 
end
dm = pi / ( u - l ) * ( 1 + tantrans( p, l, u ).^2 );