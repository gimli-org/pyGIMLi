function [x,rho,eta,j,X] = regcgls(A,b,lambda,acc,jmax,x0)
% REGCGLS - Regularized CG-based Least Squares
% [X,rho,eta,j,X] = regcgls(A,b,lambda,acc,jmax,x0)
%
% Solves implicitly the regularized normal equations
% ( A' A + lambda*I )*x = A'*b
% acc - accuracy, jmax - maximum step number, x0 - starting vector
 
% Initialization.
[m,n] = size(A); 
%IN
if nargout>4, X = []; eta=[]; rho=[]; end
%OUT
if nargin<3, lambda=0.001; end
if nargin<4, acc=1e-6; end
if nargin<5, jmax=100; end
if nargin<6, x0=zeros(n,1); end
if isempty(x0), x0=zeros(n,1); end
% Prepare for CG iteration.
x = x0;
%p = A'*b;%-lambda*x;
%leave=p'*p*(acc^2);
z = b - A*x;
p = (z'*A)'-lambda*x;
normr2 = p'*p;
acc=acc*normr2;
% Iterate.
j=0;
while(j<jmax)&(normr2>acc),
  j=j+1;  
  q = A*p;
  alpha = normr2/(q'*q+lambda*p'*p);
  x  = x + alpha*p;
  z  = z - alpha*q;
  r  = (z'*A)'-lambda*x;
  normr2old=normr2;
  normr2 = r'*r;
  beta = normr2/normr2old;
  p = r + beta*p;
  if (nargout==5), 
      X = [X x];
      rho = [rho norm(z)];
      eta = [eta norm(x)]; 
  end
end
%fprintf('%d) ',j);
if nargout<5, rho=norm(z);eta=norm(x); end
