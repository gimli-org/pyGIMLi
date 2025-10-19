function dmdp=logtransdiff(p,l,u)

% LOGTRANSDIFF - Derivative of logarithmic transformation with lower/upper bound
% m = logtransdiff(p,[lowerbound[,upperbound]])
% if upperbound is not specified or zero only lowerbound is used
% if lowerbound is not specified it is set to zero
% see also LOGTRANS, LOGTRANSINV

if nargin<2, l=0; end
if nargin<3, u=l; end
dmdp = 1. / ( p - l );
if u>l, 
    dmdp = dmdp + 1. / ( u - p ); 
end