%GAUSSM Mixture of Gaussians density estimate
%
%  W = GAUSSM(A,K,R,S,M)
%  W = A*GAUSSM([],K,R,S,M);
%
% INPUT
%   A     Dataset
%   K     Number of Gaussians to use (default: 1)
%   R,S,M Regularization parameters, 0 <= R,S <= 1, see QDC
%
% OUTPUT
%   W     Mixture of Gaussians density estimate
%
% DESCRIPTION
% Estimation of a PDF for the dataset A by a Mixture of Gaussians
% procedure. Use is made of EMCLUST(A,QDC,K). Unlabeled objects are
% neglected, unless A is entirely unlabeled or double. Then all objects
% are used. If A is a multi-class crisp labeled dataset the densities are
% estimated class by class and then weighted and combined according their
% prior probabilities. In all cases, just single density estimator W is
% computed.
%
% Note that it is necessary to set the label type of A to soft labels
% (A = LABTYPE(A,'soft') in order to use the traditional EM algorithm
% based on posterior probabilities instead of using crisp labels.
% 
% The mapping W may be applied to a new dataset B using DENSITY = B*W.
%
% SEE ALSO
% DATASETS, MAPPINGS, QDC, MOGC, EMCLUST, PLOTM, TESTC

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: gaussm.m,v 1.8 2009/02/02 21:57:39 duin Exp $

function w = gaussm(a,n,r,s,dim)

	prtrace(mfilename)
	
	if nargin < 5, dim = []; end
	if nargin < 4, s = 0; end
	if nargin < 3, r = 0; end
  if (nargin < 2)
		prwarning (2,'number of Gaussians not specified, assuming 1.');
		n = 1; 
	end

	% No arguments specified: return an empty mapping.

  mapname = 'Mixture of Gaussians';
	if (nargin < 1) | (isempty(a))
		w = mapping(mfilename,{n,r,s,dim});
		w = setname(w,mapname);
		return
	end
	
	if isa(a,'prdataset')
		labname = getname(a);
	else
		labname = '';
	end
	if ((~isdataset(a) & ~isdatafile(a)) | (getsize(a,3) ~= 1 & islabtype(a,'crisp')))
		w = mclassm(a,mapping(mfilename,n),'weight');
		w = setlabels(w,labname);
		w = setname(w,mapname);
		return
	end

	[m,k] = getsize(a);
	
	if n == 1
		[U,G] = meancov(a);
		res.mean = +U;
		res.prior= 1;

        % Regularize, if requested.
        if (s > 0) || (r > 0) 
            G = (1-r-s) * G + r * diag(diag(G)) +s*mean(diag(G))*eye(size(G,1));
        end

        % If DIM < K, extract the first DIM principal components and estimate
        % the noise outside the subspace.

        if (dim < k)
            dim = min(rank(G)-1,dim);
            [eigvec,eigval] = eig(G); eigval = diag(eigval);
            [~,ind] = sort(-eigval);

            % Estimate sigma^2 as avg. eigenvalue outside subspace.
            sigma2 = mean(eigval(ind(dim+1:end)));

            % Subspace basis: first DIM eigenvectors * sqrt(eigenvalues).
            G = eigvec(:,ind(1:dim)) * diag(eigval(ind(1:dim))) * eigvec(:,ind(1:dim))' + ...
                sigma2 * eye(k);
        end
        
		res.cov  = G;
        
		w = normal_map(res,labname,k,1);
	else
		[e,v] = emclust(a,qdc([],r,s,dim),n);	
		ncomp0 = size(v.data.mean,1);	
		iter = 0;
		while (ncomp0 ~= n & iter < 5)   % repeat until exactly n components are found
	 		[e,v1] = emclust(a,qdc([],r,s,m),n);				
 			ncomp1 = size(v1.data.mean,1);
			if ncomp1 > ncomp0
				v = v1;
				ncomp0 = ncomp1;
			end
			iter = iter + 1;
		end
		res = v.data;
		res.nlab = ones(n,1); % defines that all Gaussian components have to be
		                      % combined into a single class.
		w = prmapping('normal_map','trained',res,labname,k,1);
	end
	w = setname(w,mapname);

return
