%REDUCM Reduce to minimal space
%
%  W = REDUCM(A)
%
% Ortho-normal mapping to a space in which the dataset A exactly fits.
% This is useful for datasets with more features than objects.  For the
% objects in B = A*W holds that their dimensionality is minimum, their mean
% is zero, the covariance matrix is diagonal with decreasing variances and
% the inter-object distances are equal to those of A.
%
% For this mapping just the labeled objects in A are used, unless A is
% entirely unlabeled. In that case all objects are used.
%
% See also MAPPINGS, DATASETS, NLFISHERM, KLM, PCA

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: reducm.m,v 1.5 2010/02/08 15:29:48 duin Exp $

function W = reducm(a)

		
	if nargin < 1 | isempty(a)
		W = prmapping('reducm');
		W = setname(W,'Reduction mapping');
		return
	end

	a = testdatasize(a);
	
	% Find the subspace R of dataset 'a' (actually data matrix 'b'):
	b = +cdats(a,1);
	[m,k] = size(b);	
	[R,s,v] = prsvd(b',0);
	% Map the data:
	b = b*R;
	% Find the number of non-singular dimensions
	% (can be found from svd??)
	r = rank(b);
	if r == m, r = r-1; end
	% Order the dimensions according to the variance:	
	% Overdone if we have already the svd
	G = prcov(b);
	[F V] = preig(G); 
	[v,I] = sort(-diag(V)); 
	I = I(1:r);
	% And store the result in an affine mapping:
	R = R*F(:,I);
	W = affine(R,-mean(b*F(:,I)),a);

return
