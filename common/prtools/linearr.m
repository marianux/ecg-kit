%LINEARR Linear regression
%
%     Y = LINEARR(X,LAMBDA,N)
%
% INPUT
%   X       Dataset
%   LAMBDA  Regularization parameter (default: no regularization)
%   N       Order of polynomial (optional)
%
% OUTPUT
%   Y       Linear (or higher order) regression
%
% DESCRIPTION
% Perform a linear regression on dataset X, with regularization
% parameter LAMBDA. When N is supplied, also higher order polynomials
% are possible.
%
% SEE ALSO
% RIDGER, TESTR, PLOTR, VANDERMONDEM

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
function y = linearr(x,lambda,p)

if nargin<3
	p = 1;
end
if nargin<2
	lambda = [];
end
if nargin<1 | isempty(x)
	y = prmapping(mfilename,{lambda,p});
	y = setname(y,'Linear regression');
	return
end

if ~ismapping(lambda) %training
	[n,d] = size(x);
	X = +vandermondem(x,p);
	if isempty(lambda)
		beta = prinv(X'*X)*X'*gettargets(x);
	else
		dimp = size(X,2);
		beta = prinv(X'*X + lambda*eye(dimp))*X'*gettargets(x);
	end
	W.beta = beta;
	W.n = p;
	y = prmapping(mfilename,'trained',W,1,d,1);
	y = setname(y,'Linear regression');
else
	% evaluation
	w = getdata(lambda);
	out = vandermondem(x,w.n)*w.beta;
	y = setdat(x,out);
	
end
