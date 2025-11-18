function A=mrsfwd1dblock(K,t,zvec,wc,t2,thk)

% A=mrsfwd1dblock(K,t,zvec,wc,t2,thk)

% wc=[0.05 1 0.05];
% t2=[0.5 1.5 0.05];
% thk=[3 5];
%%
zhk=cumsum(thk);
nlay=max([length(wc) length(t2) length(thk)+1]);
nzvec=length(zvec);
%% find indices for z weight
izvec=0;
rzvec=zeros(nlay,1);
for i=1:length(thk),
    ii=find(zvec<zhk(i),1,'last');
    izvec(i+1)=ii;
    if ii<nzvec,
        rzvec(i+1)=(zhk(i)-zvec(ii))/(zvec(ii+1)-zvec(ii));
    end
end
% zvec(izvec+1)
izvec(end+1)=length(zvec)+1;
wcvec=zeros(length(zvec),1);
A=zeros(length(t),size(K,1));
for i=1:nlay,
    wcvec(:)=0;
    wcvec(izvec(i)+1:izvec(i+1)-1)=wc(i);
    if izvec(i+1)<nzvec, wcvec(izvec(i+1))=wc(i)*rzvec(i+1); end
    if izvec(i)>0, wcvec(izvec(i))=wc(i)*(1-rzvec(i)); end
    amps=K*wcvec;
    ett=exp(-t/t2(i));
    for j=1:length(amps), A(:,j)=A(:,j)+ett*amps(j); end        
end
% imagesc(log10(abs(A')))