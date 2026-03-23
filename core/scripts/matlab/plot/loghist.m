function mima=loghist(data,xli,nn)

% LOGHIST - Make histogram based on logarithmized values
% loghist(field[,bins])

if nargin<3, nn=30; end
if nargin<2, xli=[]; end
hist(log10(data),nn);
if ~isempty(xli), 
    xlim(log10([xli])); 
else
    axis tight
end
xt=rndig(10.^get(gca,'XTick'));
xtl=num2strcell(xt);
set(gca,'XTickLabel',xtl);
