function [x, w] = gauleg(x1, x2, n)
% GAULEG - Gauss-Legendre-Koefficients and Weights
% [x,w] = gauleg(x1,x2,n)
% Given the lower and upper limits of integration x1 and x2 and given n,
% this routine returns arrays x(1..n) and w(1..n) of length n,
% containing the abscissas and weights of the Gauss-Legendre
% n-point quadrature formula.
   
% For a description of the following routines see
% Numerical Recipes, Press et al.

EPS = 3.0e-6;
	
m=(n+1)/2;
xm=0.5*(x2+x1);
xl=0.5*(x2-x1);
for i=1:m
   z=cos(3.141592653*(i-0.25)/(n+0.5));   
% Starting with the above approximation to the ith root, we enter
% the main loop of refinements by Newton's method
	z1 = z + 2*EPS;
	while abs(z-z1) > EPS
      p1 = 1.0; 
      p2 = 0.0;
      for j=1:n
         p3 = p2; 
         p2 = p1;
         p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      end      
% p1 is now the desired Legendre polynomial. We next compute pp,
% its derivative, by a standard relation involving also p2, the
% polynomial of one lower order
		pp = n*(z*p1-p2)/(z*z-1.0);
		z1 = z;
      z = z1-p1/pp; % Newtons method
   end
% Scale the root to the desired interval, and put in its
% symmetric counterpart
	x(i) = xm-xl*z;
	x(n+1-i) = xm+xl*z;
% Compute the weight and ist symmetric counterpart
	w(i) = 2.0*xl/((1.0-z*z)*pp*pp);
	w(n+1-i) = w(i);
end