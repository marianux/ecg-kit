%QDC Quadratic Bayes Normal Classifier (Bayes-Normal-2)
%
%   W = QDC(A,R,S)
%
% INPUT
%   A    Dataset
%   R,S	 Regularization parameters, 0 <= R,S <= 1 
%        (optional; default: no regularization, i.e. R,S = 0)
%   M    Dimension of subspace structure in covariance matrix (default: K,
%           all dimensions)
%
% OUTPUT
%   W    Quadratic Bayes Normal Classifier mapping
%
% DESCRIPTION
% Computation of the quadratic classifier between the classes of the dataset
% A assuming normal densities. R and S (0 <= R,S <= 1) are regularization
% parameters used for finding the covariance matrix by
% 
%   G = (1-R-S)*G + R*diag(diag(G)) + S*mean(diag(G))*eye(size(G,1))
%
% This covariance matrix is then decomposed as G = W*W' + sigma^2 * eye(K),
% where W is a K x M matrix containing the M leading principal components.
% 
% The use of soft labels is supported. The classification A*W is computed by
% NORMAL_MAP.
%
% EXAMPLES
% See PREX_MCPLOT, PREX_PLOTC.
%
% REFERENCES
% 1. R.O. Duda, P.E. Hart, and D.G. Stork, Pattern classification, 2nd
% edition, John Wiley and Sons, New York, 2001. 
% 2. A. Webb, Statistical Pattern Recognition, John Wiley & Sons, 
% New York, 2002.
%
% SEE ALSO
% MAPPINGS, DATASETS, NMC, NMSC, LDC, UDC, QUADRC, NORMAL_MAP

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: qdc.m,v 1.11 2003/11/22 23:19:26 bob Exp $

function w = qdc_new(a,r,s,dim)

	prtrace(mfilename);

  if (nargin < 4)
    prwarning(4,'subspace dimensionality M not given, assuming K');
    dim = [];
  end
	if (nargin < 3)
		prwarning(4,'Regularisation parameter S not given, assuming 0.');
		s = 0; 
	end
	if (nargin < 2)
		prwarning(4,'Regularisation parameter R not given, assuming 0.');
		r = 0;
	end

	% No input arguments: return an untrained mapping.

	if (nargin < 1) || (isempty(a))
		w = mapping(mfilename,{r,s,dim});
		w = setname(w,'New_Bayes-Normal-2');
		return
	end

	islabtype(a,'crisp','soft');
	isvaldset(a,2,2); % at least 2 objects per class, 2 classes

	[m,k,c] = getsize(a);

	% If the subspace dimensionality is not given, set it to all dimensions.

	if (isempty(dim)), dim = k; end;

  if (dim < 1) || (dim > k)
    error ('Number of dimensions M should lie in the range [1,K].');
  end

	% Assert A has the right labtype.

	islabtype(a,'crisp','soft');

	[U,G] = meancov_new(a, 0);

	% Calculate means and priors.

	pars.mean  = +U;
	pars.prior = getprior(a);

	% Calculate class covariance matrices.

	pars.cov   = zeros(k,k,c);
	for j = 1:c
		F = G(:,:,j);
		
		% Regularize, if requested.

		if (s > 0) || (r > 0) 
			F = (1-r-s) * F + r * diag(diag(F)) +s*mean(diag(F))*eye(size(F,1));
		end

    % If DIM < K, extract the first DIM principal components and estimate
    % the noise outside the subspace.

        if (dim < k)
          [eigvec,eigval] = eig(F); eigval = diag(eigval);
          [dummy,ind] = sort(-eigval);

          % Estimate sigma^2 as avg. eigenvalue outside subspace.
          sigma2 = mean(eigval(ind(dim+1:end)));

          % Subspace basis: first DIM eigenvectors * sqrt(eigenvalues).
          F = eigvec(:,ind(1:dim)) * diag(eigval(ind(1:dim))) * eigvec(:,ind(1:dim))' + ...
                sigma2 * eye(k);
        end

		pars.cov(:,:,j) = F;
        
	end

	w = mapping('normal_map_new','trained',pars,getlab(U),k,c);
	w = setname(w,'Bayes-Normal-2');
	w = setcost(w,a);
    
    cFeaturesDomain = getfeatdom(a);

    w = setuser(w,cFeaturesDomain);

return;
