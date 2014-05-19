%PRRANK Call to RANK() including PRWAITBAR 
%
%	B = PRRANK(A,tol)
%
% This calls B = RANK(A,tol) and includes a message to PRWAITBAR
% in case of a large A

function B = prrank(varargin)

[m,n] = size(varargin{1});
if min([m,n]) >= 500
	prwaitbaronce('Rank of %i x %i matrix ...',[m,n]);
	B = rank(varargin{:});
	prwaitbar(0);
else
	B = rank(varargin{:});
end