function erg=logdroptrans(vals,drop)

% LOGDROPTRANS - Sign-conservative log transformation with drop tolerance

if nargin<2, 
    drop=1e-6; 
end
aa = abs( vals / drop );
aa( aa<1 ) = 1;
erg = log10( aa ) .* sign( vals );