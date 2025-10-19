function [X,rho,eta,j] = cglspar(A,b,lambda)
% REGCGLS - Regularized Least Squares with Conjugate Gradients
% [X,rho,eta,j,X] = regcgls(A,b,lambda,acc,jmax,x0)
%
% Performs k steps of the conjugate gradient algorithm applied
% implicitly to the normal equations (A'*A+lambda*I)*x = A'*b.
 
% Initialization.
if nargin<3, lambda=0.1.^(0:5); end
[m,n] = size(A); 
kmax=length(lambda);
% Prepare for CG iteration.
x=zeros(n,1);
X=zeros(n,kmax);
z = b - A*x;
p = (z'*A)';  % wegen x0=0
Z=zeros(m,kmax); % residuum of Ax=b
P=X; % update vectors
for k=1:kmax,
    Z(:,k)=z;
    P(:,k)=p-lambda(k)*x;
end
R=P; % residuum of normal equation A'Ax=A'b
normr2 = ones(kmax,1)*(p'*p);
% Iterate.
abbr=p'*p*1e-7;
j=0;
kmin=1;
wb=waitbar(0,'CG Least Squares Parallel...');
while(kmin<=kmax),
  j=j+1;  
  Q = A*P;
  normr2old=normr2;
  for k=kmin:kmax,
      q=Q(:,k);p=P(:,k);
      alpha = normr2(k)/(Q(:,k)'*Q(:,k)+lambda(k)*P(:,k)'*P(:,k));
      X(:,k)  = X(:,k) + alpha*P(:,k);
      Z(:,k)  = Z(:,k) - alpha*Q(:,k);
  end
  R(:,kmin:kmax)=(Z(:,kmin:kmax)'*A)';
  for k=kmin:kmax,
      %R(:,k)  = (Z(:,k)'*A)'-lambda(k)*X(:,k);
      R(:,k)  = R(:,k)-lambda(k)*X(:,k);
      normr2(k) = R(:,k)'*R(:,k);
      beta = normr2(k)/normr2old(k);
      P(:,k) = R(:,k) + beta*P(:,k);
  end
  okmin=kmin;
  kmin=min(find(normr2>abbr));
  if isempty(kmin), kmin=kmax+1; end
  if kmin>okmin, waitbar(kmin/kmax,wb); end
end
close(wb);
if nargout>1, 
    for k=1:kmax,
        rho(k)=norm(Z(:,k));
        eta(k)=norm(X(:,k));
    end
end
