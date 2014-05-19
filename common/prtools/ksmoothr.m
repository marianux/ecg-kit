%KSMOOTHR Kernel smoother
%
%      W = KSMOOTHR(X,H)
%
% INPUT
%   X    Regression dataset
%   H    Width parameter (default H=1)
%
% OUTPUT
%   W    Kernel smoother mapping
%
% DESCRIPTION
% Train a kernel smoothing W on data X, with width parameter H.
%
% SEE ALSO
%  KNNR, TESTR, PLOTR

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
function y = ksmoothr(x,h)

if nargin<2
	h = 1;
end
if nargin<1 | isempty(x)
	y = prmapping(mfilename,{h});
	y = setname(y,'Kernel smoother');
	return
end

if ~ismapping(h) %training: just store the training data
	[n,d] = size(x);
	W.x = +x;
	W.y = gettargets(x);
	W.h = h;
	y = prmapping(mfilename,'trained',W,1,d,1);
	y = setname(y,'Kernel smoother');
else
	% evaluation
	W = getdata(h);
	[n,d] = size(x);
	m = size(W.x,1);
	xtst = +x;
	gamma = -1/(W.h*W.h); % tiny speedup
	% now go through all test data:
	y = zeros(n,1);
	K = exp(gamma*distm(xtst,W.x));
	y = (K*W.y)./sum(K,2);
	y = setdat(x,y);
	
end
