function ss=int2hum(nn,istex)

% INT2HUM - Plot integer/double human readable

if nargin<2, istex=0; end
if nn>1,
    pre='kMGT';    
    l=0;
    while nn>1000,
        nn=nn/1000;
        l=l+1;
    end
    if nn>20, ss=sprintf('%d',round(nn)); else ss=sprintf('%.1f',nn); end
else
    pre='mµnp';
    l=0;
    while nn<1,
        nn=nn*1000;
        l=l+1;
    end
    if nn>100, ss=sprintf('%d',round(nn)); else ss=sprintf('%.1f',nn); end    
end
if l>0, ss=[ss pre(l)]; end
% if istex&&(ss(end)=='µ'), ss(end:end+2)='\mu'; end