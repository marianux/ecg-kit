%TESTR MSE for regression
%
%      E = TESTR(X,W,TYPE)
%      E = TESTR(X*W,TYPE)
%      E = X*W*TESTR([],TYPE)
%
% INPUT
%   X    Regression dataset
%   W    Regression mapping
%   TYPE Type of error measure, default: mean squared error
%
% OUTPUT
%   E    Mean squared error
%
% DESCRIPTION
% Compute the error of regression W on dataset X. The following error
% measures have been defined for TYPE:
% 'mse'    mean squared error
% 'mad'    mean absolute deviation
%
% SEE ALSO
%  RSQUARED, TESTC

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
function e = testr(x,w,type)

if nargin<3
	type = 'mse';
end
if nargin<2
	w = [];
end
if isstr(w)
	type = w;
	w = [];
end
if nargin<1 | isempty(x)
	e = prmapping(mfilename,'fixed',{w,type});
	e = setname(e,type);
	return
end

if (ismapping(w) & istrained(w))
	a = a*w;
end
switch type
case 'mse'
	e = mean((+x(:,1) - gettargets(x)).^2);
case 'mad'
	e = mean(abs(+x(:,1) - gettargets(x)));
otherwise
	error('Error %s is not implemented.',type);
end

if nargout==0
	%display results on the screen:
	fprintf('Error on %d objects: %f.\n',...
		size(x,1), e);
	clear e;
end

