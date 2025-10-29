function alfa=senscov2alpha(senscov)

% SENSCOV2ALPHA - Transform sensitivity coverage to
%                 alpha values using histogram
% alfa = senscov2alpha( senscov )

[nn,hh]=hist(senscov(:),50);
nnn=cumsum(nn)/length(senscov);
imi=min(find(nnn>0.02));
if isempty(imi), imi=1; end
ima=max(find(nnn<0.5));
if isempty(ima), ima=length(hh); end
mi=hh(imi);
ma=hh(ima);
alfa=(senscov-mi)/(ma-mi);
alfa(alfa<0)=0;
alfa(alfa>1)=1;
