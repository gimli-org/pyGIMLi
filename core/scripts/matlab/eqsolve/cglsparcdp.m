function [X,rho,eta,j] = cglsparcdp(A,b,L,D,PP,lambda)

% CGLSPARCDP - Least squares basing on conjugate gradients
%              for parallel solutions for different lambda
% solves (A'D'DA + lam*C) (P\x) = A'D'D b
%  = cglsparcdp(A,b,C,D,P,lam)
%     C .. constraint matrix
%     D .. data weighting matrix
%     P .. parameter mapping matrix
%   lam .. regularization parameters
%     X .. solution matrix (columns correspond to lams)
% [X,rho,eta] = ... yields (inexact) model and data norms
%               which can be used for L-curve criterion

% Initialization.
if nargin<5, lambda=0.1.^(0:5); end
[m,n] = size(A); 
x0=zeros(n,1);
kmax=length(lambda);
% Prepare for CG iteration.
PI=PP';
su=sum(PP,1);
we=ones(size(su));
fi=find(su);we(fi)=1./su(fi);
PI=spdiags(we(:),0,length(we),length(we))*PI;
x=PI*x0;
X=zeros(length(x),kmax);P=X;
z = D*(b - A*(PP*x)); % residuum of unregularized equation
Z=zeros(length(z),kmax);
r = (z'*D*A*PP)';%-L*(PP*x); % residuum of normal equation
p=r;
normr2 = ones(kmax,1)*(r'*r);
R=zeros(length(r),kmax);
P=zeros(length(p),kmax);
for k=1:kmax,
    Z(:,k)=z;
    P(:,k)=p; %-L*x*lambda(k)
end
acc=1e-8;
abbr=normr2(1)*acc; % goal for norm(r)^2
j=0;oldf=0;
t0=clock;
kmin=1; % all lambdas
wb=waitbar(0,'CG Least Squares Parallel CD...');
while(kmin<=kmax),
  j=j+1;  
  Q = D*(A*(PP*P));
  normr2old=normr2;
  fik=find(normr2>abbr);
%   for k=kmin:kmax,
  for l=1:length(fik),
      k=fik(l);
      q=Q(:,k);p=P(:,k);
      Pp=PP*p;
      alpha = normr2(k)/(q'*q+Pp'*(L*Pp)*lambda(k));
      X(:,k)  = X(:,k) + alpha*p;
      Z(:,k)  = Z(:,k) - alpha*q;
  end
  R(:,fik)=(Z(:,fik)'*D*A*PP)'; % common part
%   R(:,kmin:kmax)=(Z(:,kmin:kmax)'*D*A*PP)'; % common part
%   for k=kmin:kmax, % only for not converged
  for l=1:length(fik),
      k=fik(l);
      r = R(:,k);
      r  = r - PI*(L*(PP*X(:,k)))*lambda(k); % additional part
      normr2(k) = r'*r; % this is really a bit tricky!
      beta = normr2(k)/normr2old(k); 
      P(:,k) = r + beta*P(:,k);
      R(:,k) = r;
  end
  okmin=kmin; %for waitbar
  kmin=min(find(normr2>abbr)); % check for convergence
  if isempty(kmin), kmin=kmax+1; end
  %   if kmin>okmin, waitbar(kmin/kmax,wb); end
  fort=1+log10(normr2(end)/abbr)/log10(acc);
  if fort>oldf+0.05,
      waitbar(fort,wb);
      oldf=fort;
  end 
end
close(wb);
message(sprintf('Solved parallel weighted normal equations in %.1fs, %d iterations,%d lambdas',...
    etime(clock,t0),j,length(lambda)));
X=PP*X;
rho=zeros(kmax,1);
eta=zeros(kmax,1);
if nargout>1, 
    for k=1:kmax,
        rho(k)=norm(Z(:,k));
        eta(k)=sqrt(X(:,k)'*(L*X(:,k)));
    end
end
clear Q PI R Z P