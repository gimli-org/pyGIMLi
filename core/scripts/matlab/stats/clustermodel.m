function [catt,AA]=clustermodel(Mesh,att,numclust)

% CLUSTERMODEL - Cluster model on unstructured mesh
% clustered_model=clustermodel(Mesh,model)

A=zeros(Mesh.ncells,Mesh.dim);
for i=1:Mesh.ncells,
    A(i,:)=mean(Mesh.node(Mesh.cell(i,:),:));
end
for i=1:size(att,2),
    if min(att(:,i)>0),
        A(:,end+1)=log10(att(:,i));
    else
        A(:,end+1)=att(:,i);
    end
end
for i=1:size(A,2), % normalization to 0..1
    mi=min(A(:,i));ma=max(A(:,i));
    A(:,i)=(A(:,i)-mi)/(ma-mi);
end
A(:,3:end)=A(:,3:end)*2; % parameters more important than position
distA = pdist(A,'euclidean');
% distA = pdist(A,'chebychev');
% linkA = linkage(distA,'complete');
linkA = linkage(distA,'average');
[mx,my] = size(linkA);
if nargin<3,
    plot(flipud(linkA(mx-20:mx,3)),'o-');
    grid on
    snum=inputdlg('Number of clusters?');
    numclust=str2num(snum{1});
end
catt=zeros(length(att),1);
for ii=1:length(numclust),
    ic=numclust(ii);
    CM = cluster(linkA,ic);%numclust);
    if size(att,2)>1, % more than 1 given
        catt(:,ii)=CM;aa=[];
        for i=1:ic,%numclust,
            fi=find(CM==i);
            for j=1:size(att,2), aa(i,j)=median(att(fi,j)); end
        end
    else
        catt(:,ii)=att;
        for i=1:ic,%numclust,
            fi=find(CM==i);
            AA(i)=median(att(fi));
            catt(fi)=AA(i);
        end
    end
    if length(numclust)>1, AA{ii}=aa; else AA=aa; end 
end
% clf;patch2dmesh(Mesh,log10(catt));