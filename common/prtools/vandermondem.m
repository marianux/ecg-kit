%VANDERMONDEM Extend data matrix
%
%    Z = VANDERMONDEM(X,N)
%
% INPUT
%    X    Data matrix
%    N    Order of the polynomail
%
% OUTPUT
%    Z    New data matrix containing X upto order N
%
% DESCRIPTION
% Construct the Vandermonde matrix Z from the original data matrix X by
% including all orders upto N. Note that also order 0 is added:
%    Z = [ones  X  X^2  X^3 ... X^N]
% This construction allows for the trivial extension of linear methods
% to obtain polynomail regressions.
%
% SEE ALSO
%   LINEARR

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
function z = vandemondem(x,n)

if nargin<2, n=1; end
if nargin<1 | isempty(x)
	z = prmapping(mfilename,'fixed',{n});
	z = setname(z,'Vandemonde map');
	return
end

% no training, just evaluation:
dat = +x;
[m,dim] = size(dat);
I = 1:dim;
z = ones(m,n*dim+1);
for i=0:(n-1)
	z(:,(i+1)*dim+I) = dat.*z(:,i*dim+I);
end
% remove the superfluous ones:
z(:,1:dim-1) = [];
z = setdat(x,z);

return
