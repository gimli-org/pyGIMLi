function p=logtransinv(m,l,u)

% LOGTRANSINV - Inverse logarithm transformation
% p = tantransinv(m[,lowerbound[,upperbound]])
% if upperbound is not specified or zero only lowerbound is used
% if lowerbound is not specified it is set to zero
% see also LOGTRANS

if nargin<2, l=0; end
if nargin<3, u=l; end
em=exp(m);
if u>l,
    p = ( em * u + l ) ./ ( 1 + em );
else
    p = em + l;
end