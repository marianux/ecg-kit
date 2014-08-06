%RSQUARED R^2 statistic
%
%      E = RSQUARED(X,W)
%      E = RSQUARED(X*W)
%      E = X*W*RSQUARED
%
% INPUT
%    X    Regression dataset
%    W    Regression mapping
%
% OUTPUT
%    E    The R^2-statistic
%
% DESCRIPTION
% Compute the R^2 statistic of regression W on dataset X.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
%  TESTR

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
function r2 = rsquared(x,w)

if nargin==0

	r2 = prmapping(mfilename,'fixed');
	return

elseif nargin==1

	y = gettargets(x);
	yhat = +x(:,1);
	meany = mean(y);
	r2 = (sum((yhat-meany).^2))/(sum((y-meany).^2));

else

	ismapping(w);
	istrained(w);

	r2 = feval(mfilename, x*w);

end

