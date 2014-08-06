%GENREGSIN Generate sinusoidal regression data
%
%     X = GENDATSIN(N,SIGMA)
%
% INPUT
%    N      Number of objects to generate
%    SIGMA  Standard deviation of the noise
%
% OUTPUT
%    X      Regression dataset
%
% DESCRIPTION
% Generate an artificial regression dataset [X,Y] with:
%
%    y = sin(4x) + noise. 
%
% where noise is Gaussian distributed with standard deviation sigma.
%
%  X = GENDATSIN(100)
% generates 100 x,y pairs with data and default noise (sigma = 0.1).
%
%  x = (0:0.01:1)';
%  y = genregsin(x,0);
% generates the true function along the x-axis, with zero noise.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% GENDATR, GENDATSINC

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
function x = gendatsin(nrx,noise)
if nargin<2
  noise = 0.1;
end

if (length(nrx)>1)
  x = nrx;
  nrx = size(x);
  x = sin(4*x) + noise*randn(nrx);
else
  x = rand(nrx,1);
  y = sin(4*x) + noise*randn(nrx,1);
end
x = gendatr(x,y);

return

