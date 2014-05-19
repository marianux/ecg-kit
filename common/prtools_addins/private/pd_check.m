%PD_CHECK Check whether the input matrix is positive definite
%
%  OK = PD_CHECK(A)
%
% INPUT
%   A   Symmetric matrix
%
% OUTPUT
%   OK  1/0 to indicate whether A is positive definite
%
% DESCRIPTION
% Symmetric matrix A is positive definite, if the diagonal of A is 
% dominant, i.e. A(i,i) > sum_{j ~= i}  A(i,j).
% Based on Golub's version of Cholesky, with non-positivity test.
%
% REFERENCES
% G. Golub and C.V. Loan, Matrix Computations, The Johns Hopkins 
% University Press, Baltimore, 1989.
