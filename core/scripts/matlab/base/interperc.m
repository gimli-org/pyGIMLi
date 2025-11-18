function [res,res2]=interperc(att,perc)

% INTERPERC - Interpercentiles
% values = interperc(vector,percentages)

if nargin<2, perc=[5 95]; end
[NN,X]=hist(att(isfinite(att)),100);
C=cumsum(NN)/sum(NN)*100;
for i=1:length(perc),
    fi=find(C>perc(i));
    if isempty(fi),
        res(i)=max(att);
    else
        res(i)=X(min(fi));%min(att);
    end
    if perc(i)<=0, res(i)=min(att); end
    if isempty(res(i)), 
        if perc(i)<50, res(i)=min(att); else res(i)=max(att); end
    end
end
if (nargout>1)&(length(perc)>1), 
    res2=res(2);
    res=res(1);
end