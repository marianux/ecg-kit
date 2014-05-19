%LASSOR LASSO regression
%
%    W = LASSOR(X,LAMBDA)
%
% INPUT
%   X       Regression dataset
%   LAMBDA  Regularization parameter
%
% OUTPUT
%   W       LASSO regression mapping
%
% DESCRIPTION
%  The 'Least Absolute Shrinkage and Selection Operator' regression,
%  using the regularization parameter LAMBDA.
%
% SEE ALSO
%  RIDGER, LINEARR, PLOTR

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
function y = lassor(x,lambda)

if nargin<2
	lambda = 1;
end
if nargin<1 | isempty(x)
	y = prmapping(mfilename,{lambda});
	y = setname(y,'LASSO regression');
	return
end

if ~ismapping(lambda) %training
	[n,d] = size(x);
	y = gettargets(x);
	W = arrfit(+x,(y-mean(y)),lambda);
	W = [mean(y); W];
	y = prmapping(mfilename,'trained',W,1,d,1);
	y = setname(y,'LASSO regression');
else
	% evaluation
	w = getdata(lambda);
	[n,d] = size(x);
	out = [ones(n,1) +x]*w;
	y = setdat(x,out);
	
end
