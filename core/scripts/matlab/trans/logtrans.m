function m=logtrans(p,l,u)

% LOGTRANS - Logarithmic transformation with lower/upper bound
% m = logtrans(p,[lower_bound[,upper_bound]])
% if upperbound is not specified or zero only lowerbound is used
% if lowerbound is not specified it is set to zero
% see also LOGTRANSINV, LOGTRANSDIFF

if nargin<2, l=0; end
if nargin<3, u=l; end
m = log( p - l );
m( p<l ) = -1/eps;
if u>l, 
    m = m - log( u - p ); 
    m( p>u ) = 1/eps;
end
