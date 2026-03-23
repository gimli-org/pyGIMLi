function [x,rho,eta,j,X] = regcgls(A,b,P,lambda,acc,jmax,x0)
% REGCGLS - Regularized Least Squares with Conjugate Gradients
% [X,rho,eta,j,X] = regcgls(A,b,P,lambda,acc,jmax,x0)
%
% Performs k steps of the conjugate gradient algorithm applied
% implicitly to the normal equations (A'*A+lambda*I)*x = A'*b.
 
% Initialization.
[m,n] = size(A); 
%IN
if nargout>4, X = []; eta=[]; rho=[]; end
%OUT
if nargin<4, lambda=0.001; end
if nargin<5, acc=1e-6; end
if nargin<6, jmax=100; end
if nargin<7, x0=zeros(n,1); end
if isempty(x0), x0=zeros(n,1); end
% Prepare for CG iteration.
x = P\x0;
%p = A'*b;%-lambda*x;
%leave=p'*p*(acc^2);
z = b - A*(P*x);
%p = A'*z-lambda*x;
p = (z'*A*P)'-lambda*x;
normr2 = p'*p;
% Iterate.
j=0;
while(j<jmax)&(normr2>acc),
  j=j+1;  
  q = A*(P*p);
  alpha = normr2/(q'*q+lambda*p'*p);
  x  = x + alpha*p;
  z  = z - alpha*q;
  r  = (z'*A*P)'-lambda*x;
  normr2old=normr2;
  normr2 = r'*r;
  beta = normr2/normr2old;
  p = r + beta*p;
  if (nargout==5), 
      X = [X P*x];
      rho = [rho norm(z)];
      eta = [eta norm(x)]; 
  end
end
%fprintf('%d) ',j);
if nargout<5, rho=norm(z);eta=norm(x); end
x=P*x;