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

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: expandd.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function a = expandd(x,m)

	x = x(:);
	n = length(x);
	if (nargin == 1) 
		m = max(x); 
	end
	a = (x(:,ones(1,m)) == ones(n,1) * linspace(1,m,m));

return;
