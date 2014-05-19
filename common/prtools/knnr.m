%KNNR Nearest neighbor regression
%
%    Y = KNNR(X,K)
%
% INPUT
%   X    Regression dataset
%   K    number of neighbors (default K=3)
%
% OUTPUT
%   Y    k-nearest neighbor regression
%
% DESCRIPTION
% Define a k-Nearest neighbor regression on dataset X.
%
% SEE ALSO
%  LINEARR, TESTR, PLOTR

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
function y = knnr(x,k)

if nargin<2
	k = 3;
end
if nargin<1 | isempty(x)
	y = prmapping(mfilename,{k});
	y = setname(y,'k-nearest neighbor regression');
	return
end

if ~ismapping(k) %training
	[n,d] = size(x);
	W.x = +x;
	W.y = gettargets(x);
	W.k = k;
	y = prmapping(mfilename,'trained',W,1,d,1);
	y = setname(y,'k-nearest neighbor regression');
else
	% evaluation
	w = getdata(k);
	[n,d] = size(x);
	D = distm(+x,w.x);
	[sD,I] = sort(D,2);
	if n==1
		out = mean(w.y(I(:,1:w.k)));
	else
		out = mean(w.y(I(:,1:w.k)),2);
	end
	y = setdat(x,out);
	
end
