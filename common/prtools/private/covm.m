%COVM Compute covariance matrix for large datasets
% 
% 	C = COVM(A)
% 
% Similar to C = COV(A) this routine computes the covariance matrix 
% for the datavectors stored in the rows of A. No large intermediate 
% matrices are created. If class(A) is 'prdataset' then class(C) is 
% 'double'.

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: covm.m,v 1.3 2010/02/08 15:34:14 duin Exp $

function c = covm(a,n)
		[m,k] = size(a);
if nargin < 2, n = 0; end
if n ~= 1 & n ~= 0
	error('Second parameter should be either 0 or 1')
end
[loops,n0,n1] = prmem(m,k);
if loops == 1
	c = prcov(+a,n);
	c = (c+c')/2;
	return
end
c = zeros(k,k);
u = ones(n0,1)*mean(a);
for j = 1:loops
	if j == loops, n = n1; else n = n0; end
	nn = (j-1)*n0;
	b = +a(nn+1:nn+n,:) - u(1:n,:);
	c = c + b'*b;
end
c = (c + c')/2;
if n
	c = c/m;
else
	c = c/(m-1);
end
