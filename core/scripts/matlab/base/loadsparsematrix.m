function S = loadsparsematrix(sensname)

% LOADSENS - Load sparse matrix from file
% S = loadsparsematrix(sensname)

A=textread(sensname);
S=sparse(A(:,1)+1,A(:,2)+1,A(:,3));