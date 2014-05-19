%EXPANDD Expand integer vector to a matrix occurrence table
% 
%   A = EXPANDD(X,M)
% 
% INPUT
%   X   Vector containing positive integers
%   M   Dimensionality of the returned matrix A
%
% OUTPUT
%   A   Matrix 
%
% DESCRIPTION
% The vector X (e.g. numeric labels obtained from RENUMLAB) is expanded to 
% a matrix A of the size [LENGTH(X) x M] such that A(i,j) = 1 if X(i) == j 
% and A(i,j) = 0, otherwise.
% As a result SUM(A) is a frequency table of j in X. MAX(A) is a 0/1 vector 
% indicating the occurrence of some j in X. FIND(MAX(A)) gives all j that 
% occur in X.
% 
% SEE ALSO 
% RENUMLAB
