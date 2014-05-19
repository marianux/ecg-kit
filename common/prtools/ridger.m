%RIDGER Ridge Regression
%
%     W = RIDGER(X,LAMBDA)
%
% INPUT
%   X        Regression dataset
%   LAMBDA   Regularization parameter (default LAMBDA=1)
%
% OUTPUT
%   W        Ridge regression mapping
%
% DESCRIPTION
% Perform a ridge regression on dataset X, with the regularization
% parameter LAMBDA.
%
% SEE ALSO
%  LASSOR, PLOTR, LINEARR

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
function y = ridger(x,lambda)

if nargin<2
	lambda = 1;
end
if nargin<1 | isempty(x)
	y = prmapping(mfilename,{lambda});
	y = setname(y,'Ridge regression');
	return
end

if ~ismapping(lambda) %training
	[n,d] = size(x);
	y = gettargets(x);
	X = +x;
	beta = prpinv(X'*X + diag(repmat(lambda,d,1)))*X'*(y-mean(y));
	W = [mean(y); beta];   % don't forget the offset
	y = prmapping(mfilename,'trained',W,1,d,1);
	y = setname(y,'Ridge regression');
else
	% evaluation
	w = getdata(lambda);
	[n,d] = size(x);
	out = [ones(n,1) +x]*w;
	y = setdat(x,out);
	
end
