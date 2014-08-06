%GENDATR Generation of regression data
%
%      A = GENDATR(X,Y)
%
% INPUT
%   X    data matrix
%   Y    target values
%
% OUTPUT
%   A    regression dataset
%
% DESCRIPTION
% Generate a regression data from the data X and the target values Y.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
%  SCATTERR, GENDATSINC

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
function a = gendatr(x,y)

if nargin<2
	y = 1;
end
if nargin<1
	x = 0;
end
% check the sizes
[n,dim] = size(x);
if length(y)~=n
	error('Size of X and Y do not match.');
end

% store it in the dataset:
a = prdataset(x);
a = setlabtype(a,'targets',y);
if ~isa(x,'prdataset')
	fl = {};
	for i=1:dim
		fl{i} = sprintf('x_%d',i);
	end
	a = setfeatlab(a,fl);
end

return
