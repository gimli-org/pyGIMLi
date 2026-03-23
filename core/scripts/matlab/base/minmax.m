function mima=minmax(A)

% MINMAX - return minimum&maximum values of array or matrix
%          useful for quick check ranges and colorscales
% mima = minmax(A)

mima=[min(A(:)) max(A(:))];
