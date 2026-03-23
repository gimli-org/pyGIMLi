function K = kronecker(A,B)
%KRONECKER   Kronecker tensor product.
%   KRONECKER(X,Y) is the Kronecker tensor product of X and Y.
%   The result is a large matrix formed by taking all possible
%   products between the elements of X and those of Y. For
%   example, if X is 2 by 3, then KRONECKER(X,Y) is
%
%      [ X(1,1)*Y  X(1,2)*Y  X(1,3)*Y
%        X(2,1)*Y  X(2,2)*Y  X(2,3)*Y ]
%
%   If either X or Y is sparse, only nonzero elements are multiplied
%   in the computation, and the result is sparse.
%
%   Class support for inputs X,Y:
%      float: double, single
%
%   NOTE: This function does exactly what Matlab KRON does, but for large
%      full matrices, the engine uses BSXFUN to accelerate the calculation.   
%      Another advantage is no intermediates large matrices is generated
%      (four temporary arrays in case of KRON)
%
%   Benchmark on Intel Core2 Duo T7250 @2GHz and 2Go RAM
%   Size A/B  Speed gain
%      10       1.17    
%      20       3.48    
%      30       3.78    
%      40       3.73    
%      50       3.68    
%      60       4.22    
%      70       3.81
%
%   Restriction: MATLAB 2007A or later is required
%
%   See also: KRON
%
%   Author: Bruno Luong <brunoluong@yahoo.com>
%   History:
%       Original 21-Jun-2009



if ~issparse(A) && ~issparse(B)
    if ndims(A) > 2 || ndims(B) > 2
        error('kronecker:TwoDInput','Inputs must be 2-D.');
    end
    % Both inputs are full, result is full. This is faster than
    % MATLAB stock kron (indexing based)
    [ma na] = size(A);
    [mb nb] = size(B);
    A = reshape(A,[1 ma 1 na]);
    B = reshape(B,[mb 1 nb 1]);
    K = bsxfun(@times,A,B);
    K = reshape(K,[ma*mb na*nb]);
    
else % One of the input matrix is sparse
    
    % Call MATLAB stock KRON
    K = kron(A,B);
    
end
