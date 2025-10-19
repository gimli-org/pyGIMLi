function [X,flag,relres,iter] = cholpcg(A,B,tol,maxit,M1,M2,x0)

% CHOLPCG - Conjugate Gradient method with incomplete Cholesky factorization
% [x,flag] = cholpcg(A,b,tol,maxit,CholA,x0);

if (nargin < 2)
   error('Not enough input arguments.');
end

% Assign default values to unspecified parameters
if (nargin < 3) | isempty(tol)
   tol = 1e-6;
end
if (nargin < 4) | isempty(maxit)
   maxit = min(n,20);
end
rhs=size(B,2);
n=length(A);
% Check for all zero right hand side vector => all zero solution
n2b = sqrt(sum(B.^2))';                      % Norm of rhs vector, b
X = zeros(length(A),rhs);                  % then  solution is all zeros
if (max(n2b) == 0)                       % if    rhs vector is all zeros
   flag = 0;                        % a valid solution has been obtained
   relres = 0;                      % the relative residual is actually 0/0
   iter = 0;                        % no iterations need be performed
   return
end

% Set up for the method
flag = 1;
Xmin = X;                          % Iterate which has minimal residual so far
imin = 0;                          % Iteration at which xmin was computed
tolb = tol * n2b;                  % Relative tolerance
R = B - A * X;                  % Zero-th residual
normr = sqrt(sum(R.^2))';                   % Norm of residual

if (normr <= tolb)                 % Initial guess is a good enough solution
   flag = 0;
   relres = Normr ./ n2b;
   iter = 0;
   return
end

normrmin = normr;                  % Norm of minimum residual
rho = zeros(rhs,1);
stag = 0;                          % stagnation of the method

% loop over maxit iterations (unless convergence or failure)

for i = 1 : maxit
    i
   Y = M1 \ R;
   if isinf(max(abs(Y(:))))
      flag = 2;
      break
   end
   
   Z = M2 \ R;
   if isinf(max(abs(Z(:))))
      flag = 2;
      break
   end
   
   rho1 = rho;
   rho = sum(R.*Z)';
   if find((rho == 0) | isinf(rho))
      flag = 4;
      break
   end
   if (i == 1)
      P = Z;
   else
      beta = rho ./ rho1;
      if find((beta == 0) | isinf(beta))
         flag = 4;
         break
      end
      for k=1:rhs
        P(:,k) = Z(:,k) + P(:,k)*beta(k);
    end
   end
   Q = A * P;
   pq = sum(P.*Q)'; %p' * q;
   if find((pq <= 0) | isinf(pq))
      flag = 4;
      break
   else
      alpha = rho ./ pq;
   end
   if isinf(alpha)
      flag = 4;
      break
   end
   if (alpha == 0)                  % stagnation of the method
      stag = 1;
   end
   
   % Check for stagnation of the method !!! still to do
%    if (stag == 0)
%       stagtest = zeros(n,1);
%       ind = (x ~= 0);
%       stagtest(ind) = p(ind) ./ x(ind);
%       stagtest(~ind & p ~= 0) = Inf;
%       if (abs(alpha)*norm(stagtest,inf) < eps)
%          stag = 1;
%       end
%    end
   
   for k=1:rhs,
       X(:,k) = X(:,k) + P(:,k)*alpha(k);               % form new iterate
   end
   %normr = sqrt(sum((B - A * X)^.2))';
   normr=sqrt(sum(R.^2))';
   
   if ~find(normr > tolb)               % check for convergence
      flag = 0;
      iter = i;
      break
   end
   
   if (stag == 1)
      flag = 3;
      break
   end
   
   for k=1:rhs,
       R(:,k) = R(:,k) - Q(:,k)*alpha(k);
   end
   
end                                % for i = 1 : maxit

