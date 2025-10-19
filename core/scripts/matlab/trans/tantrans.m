function m=tantrans(p,l,u)

% TANTRANS - Tangens transformation with lower and upper bound
% m = tantrans(p,lower_bound,upperbound)
% both lower and upper bound are required
% See also TANTRANSINV, TANTRANSDIFF

if nargin<3, 
    error('Specify model, lower, and upper bound!'); 
end
% m=tan(((p-l)/(u-l)-0.5)*pi);
m = - cot( ( p - l ) ./ ( u - l ) * pi );
m( p<l ) = -1/eps;
m( p>u ) = 1/eps;