function result = fcmcluster(data,param)

% FCMCLUSTER - Fuzzy c-means clustering
% result = fcmcluster(data,param)
% data or data.X (Nxn)
% param is a structure of
%      .c number of clusters
%      .m fuzziness (>=1, default=2)
%      .e tolerance (0.001)
%      .ro cluster volumes
%      .norm normalize values (0 or 1)
%      .log take logarithms for ith column if log(i)>0
%      .range limit range of values
% result is a structure of
%      .cl associated cluster number
%      .mf normalized membership function (0 to 1)
%      .data.f partition matrix (Nxc) 
%      .data.d distance matrix
%      .cluster.v cluster centers (cxn)
%      .P covariance matrix
%      .cost cost function over iterations

%data normalization
if isfield(data,'X'), X=data.X; else X=data; end
if ~isfield(param,'c'), error('specify number of clusters!'); end
f0=param.c;
%checking the parameters given
%default parameters
if isfield(param,'m'), m = param.m; else m = 2; end;
if isfield(param,'e'), e = param.m; else e = 1e-4; end;
if isfield(param,'log'),
    for i=1:size(X,2),
        if param.log(min(length(param.log),i)), 
            param.log(i)=1;
            X(:,i)=log10(X(:,i));
        end
    end
end
if isfield(param,'range')&&(size(param.range,2)==size(data,2))&&(size(param.range,1)>1),
    display('ranging values');
    
end
if isfield(param,'norm')&&(param.norm>0),
    display('normalizing values');
    for i=1:size(X,2),
        mi(i)=min(X(:,i));ma(i)=max(X(:,i));
        X(:,i)=(X(:,i)-mi(i))/(ma(i)-mi(i));
    end
end

[N,n] = size(X);
[Nf0,nf0] = size(f0); 
X1 = ones(N,1);

% Initialize fuzzy partition matrix
rand('state',0)
if max(Nf0,nf0) == 1, 		% only number of cluster given
  c = f0;
  mm = mean(X);             %mean of the data (1,n)
  aa = max(abs(X - ones(N,1)*mm)); %
  v = 2*(ones(c,1)*aa).*(rand(c,n)-0.5) + ones(c,1)*mm;
  for j = 1 : c,
    xv = X - X1*v(j,:);
    d(:,j) = sum((xv*eye(n).*xv),2);
  end;
  d = (d+1e-10).^(-1/(m-1));
  f0 = (d ./ (sum(d,2)*ones(1,c)));
  
else
  c = size(f0,2);
  fm = f0.^m; sumf = sum(fm);
  v = (fm'*X)./(sumf'*ones(1,n)); %
end;

f = zeros(N,c);                % partition matrix
iter = 0;                       % iteration counter

% Iterate
while  max(max(f0-f)) > e
  iter = iter + 1;
  f = f0;
  % Calculate centers
  fm = f.^m;
  sumf = sum(fm);
  v = (fm'*X)./(sumf'*ones(1,n));
  for j = 1 : c,
    xv = X - X1*v(j,:);
    d(:,j) = sum((xv*eye(n).*xv),2);
  end;
  distout=sqrt(d);
  J(iter) = sum(sum(f0.*d));
  % Update f0
  d = (d+1e-10).^(-1/(m-1));
  f0 = (d ./ (sum(d,2)*ones(1,c)));
end

fm = f.^m; 
sumf = sum(fm);

%results
result.data.f=f0;
result.data.d=distout;
result.cluster.v=v;
result.iter = iter;
result.cost = J;
%
[maxmf,result.cl]=max(result.data.f,[],2);
result.mf=(maxmf*param.c-1)/(param.c-1);
if isfield(param,'norm')&&(param.norm>0),
   for i=1:size(result.cluster.v,2),
       result.cluster.v(:,i)=result.cluster.v(:,i)*(ma(i)-mi(i))+mi(i);
   end
end
if isfield(param,'log'),
   for i=1:size(result.cluster.v,2),
       if param.log(i), result.cluster.v(:,i)=10.^result.cluster.v(:,i); end
   end
end

% [result.mf,result.cl]=max(result.data.f,[],2);
% result.cdef=(result.mf*param.c-1)/(param.c-1); %0-1 ranged value of membership
result.param=param;