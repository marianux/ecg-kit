%NLFISHERM Non-linear Fisher Mapping according to Marco Loog
% 
%   W = NLFISHERM(A,N)
% 
% INPUT
%   A   Dataset
%   N   Number of dimensions (optional; default: MIN(K,C)-1, where
%       K is the dimensionality of A and C is the number of classes)
%
% OUTPUT
%   W   Non-linear Fisher mapping
%
% DESCRIPTION  
% Finds a mapping of the labeled dataset A to a N-dimensional linear 
% subspace emphasizing the class separability for neighboring classes.
% 
% REFERENCES
% 1. R. Duin, M. Loog and R. Haeb-Umbach, Multi-Class Linear Feature 
% Extraction by Nonlinear PCAM, in: ICPR15, 15th Int. Conf. on Pattern 
% Recognition, vol.2, IEEE Computer Society Press, 2000, 398-401.
% 2. M. Loog, R.P.W. Duin and R. Haeb-Umbach, Multiclass Linear Dimension
% Reduction by Weighted Pairwise Fisher Criteria, IEEE Trans. on
% Pattern Analysis and Machine Intelligence, vol.23, no.7, 2001, 762-766.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, FISHERM, KLM, PCA

% Copyright: M. Loog, R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: nlfisherm.m,v 1.4 2010/02/08 15:29:48 duin Exp $

function W = nlfisherm(a,n)
		if (nargin < 2)  
		n = []; 
	end

	% No input data, an untrained mapping returned.
	if (nargin < 1) | (isempty(a))
		W = prmapping('nlfisherm',n);
		W = setname(W,'Non-linear Fisher mapping');
		return;
	end

	islabtype(a,'crisp');
	isvaldfile(a,1,2); % at least 2 objects per class, 2 classes
	a = testdatasize(a);
	
	[m,k,c] = getsize(a);
	prior = getprior(a);
	a = setprior(a,prior);
	if (isempty(n))
		n = min(k,c)-1;
		prwarning(4,'Dimensionality N not supplied, assuming MIN(K,C)-1.');
	end

	if (n >= m) | (n >= c)
		error('Dataset too small or has too few classes for demanded output dimensionality.')
	end

	% Non-linear Fisher mapping is determined by the eigenvectors of CW^{-1}*CB,
  % where CW is the within-scatter, understood as the averaged covariance 
	% matrix weighted by the prior probabilities, and CB is the between-scatter,
	% modified in a nonlinear way.
	% To simplify the computations, CW can be set to the identity matrix.
	w = klms(a);						 
	% A is changed such that CW = I and the mean of A is shifted to the origin.
	b = a*w;							
	k = size(b,2);
	u = +meancov(b);
	d = +distm(u);				% D is the Mahalanobis distance between the classes.  

	% Compute the weights E to be used in the modified between-scatter matrix G
	% E should diminish the influence of large distances D.
	e = 0.5*erf(sqrt(d)/(2*sqrt(2))); 
	G = zeros(k,k);
	for j = 1:c
		for i=j+1:c
			% Marco-Loog Mapping
			G = G + prior(i)*prior(j)*e(i,j)*(u(j,:)-u(i,:))'*(u(j,:)-u(i,:))/d(i,j); 
			G = (G + G')/2;		% Avoid a numerical inaccuracy: cov. matrix should be symmetric!
		end
	end

	% Perform the eigendecomposition of the modified between-scatter matrix.
	[F,V] = preig(G); 			
	[v,I] = sort(-diag(V)); 
	I = I(1:n);
	rot = F(:,I);
	off = -mean(b*F(:,I));

	% After non-linear transformations, NLFISHERM is stored as an affine (linear) map.
	W = affine(rot,off,a);
	W = setname(w*W,'Non-linear Fisher mapping');

return;
