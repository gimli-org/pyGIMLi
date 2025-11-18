function p=tantransinv(m,l,u)

% TANTRANSINV - Inverse tangens transformation with lower and upper bound
% p = tantransinv(m,lowerbound,upperbound)
% both lower and upper bound are required
% See also TANTRANS

if nargin<3, 
    error('Specify model, lower, and upper bound!'); 
end
p = atan(m) / pi * (u-l) + ( l + u ) / 2;