%PRCOV Call to COV() including PRWAITBAR 
%
%	B = PRCOV(A,tol)
%
% This calls C = COV(A, ...) and includes a message to PRWAITBAR
% in case of a large A

function B = prcov(varargin)

[m,n] = size(varargin{1});
if m*n*n > 1e9
	prwaitbaronce('covariance of %i x %i matrix ...',[m,n]);
	B = cov(varargin{:});
	prwaitbar(0);
else
	B = cov(varargin{:});
end