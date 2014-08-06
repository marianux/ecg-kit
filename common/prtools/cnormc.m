%CNORMC Classifier normalisation for ML posterior probabilities
% 
%   W = CNORMC(W,A)
% 
% INPUT
% 	W  Classifier mapping
%		A  Labeled dataset
%
% OUTPUT
%		W  Scaled classifier mapping
%
% DESCRIPTION
% The mapping W is scaled such that the likelihood of the posterior
% probabilities of the samples in A, estimated by A*W*SIGM, are maximized.
% This is particularly suitable for two-class discriminants. To obtain
% consistent classifiers in PRTools, it is necessary to call CNORMC in the
% construction of all classifiers that output distances instead of densities
% or posterior probability estimates.
%
% If A has soft labels or target labels, W is returned without change.
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, CLASSC

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: cnormc.m,v 1.5 2008/07/03 09:14:02 duin Exp $

function w = cnormc (w,a)

		
	% Parameters
	reg     = 1e-7; 		% Regularization, important for non-overlapping classes.
	epsilon = 1e-6;			% Stop when change in likelihood falls below this.

	% Check arguments.
	if (~ismapping(w))
		error('first argument should be a mapping');
	end
	if (islabtype(a,'targets'))
		prwarning(4,'Dataset has target labels, not taking any action.');
		return;
	end
	if (islabtype(a,'soft'))
		prwarning(4,'Dataset has soft labels, scaling optimization is not yet implemented.');
		w = setscale(w,0.01); % just a trial
		w = setout_conv(w,1);
		return;
	end

	[m,k,c] = getsize(a); nlab = getnlab(a); p = getprior(a,0);
	mapped = double(a*w);						% Calculate mapping outputs.
	if (c == 2) & (size(mapped,2) == 1)		
		mapped = [mapped -mapped];				% Old (perhaps unnecessary) transformation 
														% for classifiers that give distances.
	end
	mapped = mapped(m*(nlab-1)+[1:m]');	% Select output corresponding to label
																			% (uses 1D matrix addressing).

	% Initialise variables.
	scale = 1e-10; likelihood = -inf; likelihood_new = -realmax;

	% Iteratively update scale to maximize likelihood of outputs.
	while (abs(likelihood_new - likelihood) > epsilon)
		
		% Class A is, for each object, the class the object is classified as;
		% class B stands for all other classes.
		pax = sigm(mapped*scale,1);  			% Estimate of P(class A|x).
		pbx = 1 - pax; 							% Estimate of P(class B|x).
		likelihood = likelihood_new; 
		likelihood_new = mean(p(nlab)'.*log(pax + realmin)) - reg*log(scale+1);

		% Gradient step to maximize likelihood.
		scale = scale + ((p(nlab)'.*pbx)' * mapped - reg*m/(scale+1)) / ...
								    ((p(nlab)'.*mapped.*pax)'*(mapped.*pbx) - ...
											m*reg/((scale+1)*(scale+1)) + realmin);
	end

	% Dirty trick to avoid 0 and inf.
	scale = min(max(scale,0.1/(std(mapped)+1e-16)),1e16);

	% Adjust mapping: set scale and set output conversion to SIGM(1).
	w = setscale(w,scale);
	w = setout_conv(w,1);

return;

