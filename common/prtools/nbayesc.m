%NBAYESC Bayes Classifier for given normal densities
% 
%   W = NBAYESC(U,G)
% 
% INPUT
%   U  Dataset of means of classes   
%   G  Covariance matrices (optional; default: identity matrices)
%
% OUTPUT
% 	W  Bayes classifier
%
% DESCRIPTION
% Computation of the Bayes normal classifier between a set of classes.
% The means, labels and priors are defined by the dataset U of the size
% [C x K]. The covariance matrices are stored in a matrix G of the 
% size [K x K x C], where K and C correspond to the dimensionality and 
% the number of classes, respectively. 
% 
% If C is 1, then G is treated as the common covariance matrix, yielding
% a linear solution. For G = I, the nearest mean solution is obtained.
% 
% This routine gives the exact solution for the given parameters, while
% the trainable classifiers QDC and LDC give approximate solutions, based
% on the parameter estimates from a training set. For a given dataset, U 
% and G can be computed by MEANCOV.
%
% EXAMPLES
% [U,G] = MEANCOV(GENDATB(25));
% W = NBAYESC(U,G);
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, QDC, LDC, NMC.

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: nbayesc.m,v 1.4 2009/01/04 21:11:07 duin Exp $

function W = nbayesc(U,G);

		[cu,ku] = size(U);		% CU is the number of classes and KU - the dimensionality
	if nargin == 1,
		prwarning(4,'Covariance matrix is not specified, the identity matrix is assumed.');  
		G = eye(ku);
	end

	[k1,k2,c] = size(G);	% C = 1, if G is the common covariance matrix.

	if (c ~= 1 & c ~= cu) | (k1 ~= k2) | (k1 ~= ku)
		error('Covariance matrix or a set of means has a wrong size.')
	end

	pars.mean  = +U;
	pars.cov   = G;
	pars.prior = getprior(U);

	%W = prmapping('normal_map','trained',pars,getlablist(U),ku,cu);
	W = normal_map(pars,getlablist(U),ku,cu);
	W = setname(W,'BayesNormal');
	W = setcost(W,U);

return;