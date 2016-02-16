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

% $Id: pd_check.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function ok = pd_check(a)

	% Check now added.
	if ~issym(a),
		error('A is not symmetric.');
	end

	n   = size(a,1); 
	ok  = 0; 
	tol = 1e-15;
	for j = 1:n
    if (j > 1)
			a(j:n,j) = a(j:n,j) - a(j:n,1:j-1) * a(j,1:j-1)';
		end
    if (a(j,j) < tol) 
			return;
		end
    a(j:n,j) = a(j:n,j) / sqrt(a(j,j));
	end
	ok = 1;   

%EP??? I am not sure how to understand the lines below.
	% ok = 1 if a is safely pos.def, ie each diag > tol in chol.faktor
	% otherwise ok = 0;
	%G = tril(a);  This could be used as a Cholesky factor.

return;
