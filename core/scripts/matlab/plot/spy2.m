function spy2(A,col)
%SPY2 An improvement of SPY
%   SPY2(A)  where A is a matrix, will plot
%   the matrix with color codes.
%   Elements containing NaN or Inf will be counted as
%   non-zero but will be displayed
%   as zero elements.
%
%   SPY2(A,COLOR)  will plot zero elements with
%   white color and nonzero elements with black
%   color.
%   Elements containing NaN or Inf will be counted as
%   non-zero and will also be displayed
%   as non-zero elements.
%
%   The displayed matrix will have sides proportional to
%   the number of elements in each row/column.
%   If the matrix contains complex numbers, only the real
%   part will be displayed.
%
%   See also SPY.

% Copyright 2000-11-14, B. Rasmus Anthin.
% Revision 2003-09-15, 2003-09-16.

error(nargchk(1,2,nargin))

A=real(double(A));
A(~isfinite(A))=realmin;
A(end+1,end+1)=0;
x=1:size(A,2);
y=1:size(A,1);
[X,Y]=meshgrid(x-.5,y-.5);   %i use surf and meshgrid instead of pcolor so that
if nargin==1
   surf(X,Y,A)                  % the element numbers will be correctly aligned
   colormap gray                %just change this afterwards to get another colormap
   colorbar
else
   A=~~A;
   A(1:end,end)=0;
   A(end,1:end)=0;
   surf(X,Y,A)
   if ischar(col)
      switch(lower(col))
      case {'y','yellow'}
         col=[1 1 0];
      case {'m','magenta'}
         col=[1 0 1];
      case {'c','cyan'}
         col=[0 1 1];
      case {'r','red'}
         col=[1 0 0];
      case {'g','green'}
         col=[0 1 0];
      case {'b','blue'}
         col=[0 0 1];
      case {'w','white'}      %very bad idea :-)
         col=[1 1 1];
      case {'k','black'}
         col=[0 0 0];
      otherwise
         col=[0 0 1];         %default to blue color
      end
   end
   if isempty(col), col=[0 0 1];end
   colormap([1 1 1;col])
end
view(2)
axis ij equal tight
box on
xlabel(['nz = ' int2str(nnz(A))])