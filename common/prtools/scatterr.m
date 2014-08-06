%SCATTERR Scatter regression data
%
%      H = SCATTERR(X,CLRS)
%
% INPUT
%    X      Regression dataset
%    CLRS   Plot string (default CLRS = 'k.')
%
% OUTPUT
%    H      Vector of handles
%
% DESCRIPTION
% Scatter the regression dataset X with marker colors CLRS.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
%  PLOTR

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
function h = plottarget(x,clrs)
if nargin<2
	clrs = 'k.';
end

y = gettargets(x);
h = plot(+x(:,1),y,clrs);
fl = getfeatlab(x);
xlabel(fl(1,:));
ylabel('target');
% also set the identifiers:
ud = get(h,'UserData');
ud.ident = getident(x);
set(h,'UserData',ud);

if nargout<1
	clear h;
end

return
