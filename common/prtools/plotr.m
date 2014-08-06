%PLOTR Plot regression
%
%    PLOTR(W)
%    PLOTR(W,CLR)
%
% Plot the regression function W, optionally using plot string CLR.
% This plot string can be anything that is defined in plot.m.
% For the best results (concerning the definition of the axis for
% instance) it is wise to first scatter the regression data using
% SCATTERR(A).
%
% The resolution of the plot is determined by the global parameter
% GRIDSIZE.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
%  SCATTERR, PLOTC

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
function h = plotr(w,clrs)
if nargin<2
	clrs = 'k-';
end

if ~isa(w,'prmapping') 
	error('Mapping expected. Old PLOTR has been renamed in PLOTE.');
end

% define the input on a grid:
V = axis;
g = gridsize;
%because we are plotting 1D, we can affort many more points compared to
%2D:
g = g*g;
x = linspace(V(1),V(2),g)';

% find the regression outputs:
y = x*w;

% and plot:
hold on; h = plot(x,+y,clrs);

% avoid unnecessary output:
if nargout<1, clear h; end

return
