%PRSVD Call to SVD() including PRWAITBAR 
%
%	VARARGOUT = PRSVD(VARARGIN)
%
% This calls B = RANK(A,tol) and includes a message to PRWAITBAR
% in case of a large A

function varargout = prsvd(varargin)

[m,n] = size(varargin{1});
varargout = cell(1,nargout);
if min([m,n]) >= 500
	prwaitbaronce('SVD of %i x %i matrix ...',[m,n]);
	[varargout{:}] = svd(varargin{:});
	prwaitbar(0);
else
	[varargout{:}] = svd(varargin{:});
end