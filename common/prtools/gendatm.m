%GENDATM Generation of multi-class 2-D data
% 
% 	A = GENDATM(N)
% 
% INPUT
%   N   Vector of class sizes (default: 20)
%
% OUTPUT
%   A   Dataset
%
% DESCRIPTION
% Generation of N samples in 8 classes of 2 dimensionally distributed data
% vectors. Classes have equal prior probabilities. If N is a vector of
% sizes, exactly N(I) objects are generated for class I, I = 1..8.
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, PRDATASETS

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: gendatm.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function a = gendatm(n)

	  if (nargin == 0)
		prwarning(3,'number of samples to generate not specified, assuming 20');
		n = repmat(20,1,8); 
	end;

	% Set equal priors and generate a class distribution according to it.

	p = repmat(1/8,1,8); n = genclass(n,p);

	% Generate 8 classes...

	a1 = +gendath(n(1:2));			% ...first 2 classes: Highleyman data.
	a2 = +gendatc(n(3:4))./5;		% ...next 2 classes : spherical classes.
	a3 = +gendatb(n(5:6))./5;		% ...next 2 classes : banana data.
	a4 = +gendatl(n(7:8))./5;		% ...next 2 classes : Lithuanian data.

	% Glue classes together with some proper offsets.

	a = [a1; a2+5; a3+repmat([5,0],n(5)+n(6),1); a4+repmat([0 5],n(7)+n(8),1)];

	lab = genlab(n,['a';'b';'c';'d';'e';'f';'g';'h']);
	a = prdataset(a,lab,'name','Multi-Class Problem');
  a = setprior(a,0); % make all classes equally probable

return
