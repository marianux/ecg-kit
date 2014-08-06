%GENDATSINC Generate Sinc data
%
%      A = GENDATSINC(N,SIGMA)
%
% INPUT
%    N      Number of objects to generate
%    SIGMA  Standard deviation of the noise (default SIGMA=0.1)
%
% OUTPUT
%    A      Regression dataset
%
% DESCRIPTION
%
% Generate the standard 1D Sinc data containing N objects, with Gaussian
% noise with standard deviation SIGMA. 
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
%  GENDATR, GENDATLIN, GENDATSIN

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
function a = gendatsinc(n,sig)

if nargin<2
	sig = 0.1;
end
if nargin<1
	n = 25;
end
% input data between -5 and +5
x = -5+10*rand(n,1);
% avoid problems with x==0: for x==0 the result is 1 anyway:
y = ones(size(x));
I = find(x); % find the x's unequal to 0
y(I) = sin(pi*x(I))./(pi*x(I));

a = prdataset(x);
a = setlabtype(a,'targets',y);
a = setfeatlab(a,'x_1');

return
