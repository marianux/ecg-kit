%GENDATL Generation of Lithuanian classes
% 
%  A = GENDATL(N,S)
%	
% INPUT
%  N  Number of objects per class (optional; default: [50 50])
%  S  Standard deviation for the data generation (optional; default: 1)
%
% OUTPUT
%  A  Dataset
%
% DESCRIPTION 
% Generation of Lithuanian classes, a 2-dimensional, 2-class dataset A
% of N objects according to the definition given by Raudys. 
% The data is uniformly distributed along two sausages and is superimposed
% by a normal distribution with standard deviation S in all directions. 
% Class priors are P(1) = P(2) = 0.5.
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, PRDATASETS

% Copyright: M. Skurichina, R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: gendatl.m,v 1.3 2007/06/19 11:44:14 duin Exp $

function a = gendatl(N,s)
		if nargin < 1, 
		prwarning(3,'Class cardinalities are not specified, assuming [50 50].');
		N = [50 50]; 
	end
	if nargin < 2, 
		prwarning(4,'Standard deviation for the data generation is not specified, assuming 1.');
		s = 1; 
	end
	if (length(N) == 1),
		N(2) = N(1);
	end;

	u1	= 2*pi/3*(rand(N(1),1)-0.5*ones(N(1),1));
	u2	= 2*pi/3*(rand(N(2),1)-0.5*ones(N(2),1));
	a 	= [[10*cos(u1) + s*randn(N(1),1) 10*sin(u1) + s*randn(N(1),1)]; ...
				[6.2*cos(u2) + s*randn(N(2),1) 6.2*sin(u2) + s*randn(N(2),1)]];
	lab = genlab(N);
	a 	= prdataset(a,lab,'name','Lithuanian Classes');
	a   = setlablist(a,[1 2]');
	a   = setprior(a,[0.5 0.5]);
return;
