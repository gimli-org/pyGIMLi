function [cc,cl,mf] = clusterplot(para1,para2,ncluster,names,marker)

% CLUSTERPLOT - two variables crossplot using cluster colors
% clusterplot(var1,var2,ncluster[,names])
% clusterplot(var1,var2,clusterresult[,names])

if nargin<5, marker='.'; end
if (nargin<4)||isempty(names), names={'resistivity','velocity'}; end
if nargin<3, ncluster=3; end
% cols={'b','br','bgr','bcgr','bcgrm','bcgyrm'};
cols={'b','br','bgr','bcyg','bcgyr','bcgyrm','bcgyrmk'};
if isstruct(ncluster), % already clustering
    result=ncluster;
else % do the clustering
    param.c=ncluster;
    result=fcmcluster([log10(para1(:)) log10(para2(:))],param);
end
%cc=10.^result.cluster.v;
col=cols{min(result.param.c,length(cols))};
ncols=length(col);
for i=1:result.param.c,
    fi=find(result.cl==i);
    coli=col(mod(i-1,ncols)+1);
    loglog(para1(fi),para2(fi),[coli marker],'MarkerSize',5);
    hold on
end
xlim(minmax(para1));
ylim(minmax(para2));
if isstruct(result)&&isfield(result,'param')&&isfield(result.param,'log'),
    if result.param.log(1)==0, set(gca,'XScale','linear'); end
    if result.param.log(min(2,length(result.param.log)))==0, set(gca,'YScale','linear'); end
end
hold off
grid on
xlabel(names{1});
ylabel(names{2});