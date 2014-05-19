%ISSYM Checks whether a matrix is symmetric
%
%   OK = ISSYM(A,DELTA)
% 
% INPUT
% 	A      Dataset
% 	DELTA  Parameter for the precision check (optional; default: 1e-12)
%
% OUTPUT
% 	OK     1 if the matrix A is symmetric and 0, otherwise.
%
% DESCRIPTION 
% A is considered as a symmetric matrix, when it is square and
% max(max(A-A')) is smaller than DELTA.
%

%
% Robert P.W. Duin, Elzbieta Pekalska, ela@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
%

function [ok,nn] = issym(A,delta)

	if nargin < 2, 
		prwarning(6,'The precision is not provided, set up to 1e-12.');
		delta = 1e-12; 
	end

	A = +A;
	[m,k] = size(A);
	if m ~= k,
	  error ('Matrix should be square.')
	end
	nn = max(max((A-A')));
	ok = (nn < delta);
return;

