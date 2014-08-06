%QDC Quadratic Bayes Normal Classifier (Bayes-Normal-2)
%
%   [W,R,S,M] = QDC(A,R,S,M)
%   [W,R,S,M]  = A*QDC([],R,S,M)
%   [W,R,S,M]  = A*QDC(R,S,M)
%
% INPUT
%   A    Dataset
%   R,S	 Regularization parameters, 0 <= R,S <= 1 
%        (optional; default: no regularization, i.e. R,S = 0)
%   M    Dimension of subspace structure in covariance matrix (default: K,
%        all dimensions)
%
% OUTPUT
%   W    Quadratic Bayes Normal Classifier mapping
%   R    Value of regularization parameter R as used 
%   S    Value of regularization parameter S as used
%   M    Value of regularization parameter M as used
%
% DESCRIPTION
% Computation of the quadratic classifier between the classes of the dataset
% A assuming normal densities. R and S (0 <= R,S <= 1) are regularization
% parameters used for finding the covariance matrix by
% 
%   G = (1-R-S)*G + R*diag(diag(G)) + S*mean(diag(G))*eye(size(G,1))
%
% This covariance matrix is then decomposed as G = W*W' + sigma^2 * eye(K),
% where W is a K x M matrix containing the M leading principal components
% and sigma^2 is the mean of the K-M smallest eigenvalues.
% 
% The use of soft labels is supported. The classification A*W is computed by
% NORMAL_MAP.
%
% If R, S or M is NaN the regularisation parameter is optimised by REGOPTC.
% The best result are usually obtained by R = 0, S = NaN, M = [], or by
% R = 0, S = 0, M = NaN (which is for problems of moderate or low dimensionality
% faster). If no regularisation is supplied a pseudo-inverse of the
% covariance matrix is used in case it is close to singular.
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
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, REGOPTC, NMC, NMSC, LDC, UDC, QUADRC, NORMAL_MAP

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: qdc.m,v 1.8 2010/02/08 15:29:48 duin Exp $

function [w,r,s,dim] = qdc(varargin)

  mapname = 'Bayes-Normal-2';
  argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],0,0,[]);
  
  if mapping_task(argin,'definition')
    
    w = define_mapping(argin,'untrained',mapname);
    
	elseif mapping_task(argin,'training')			% Train a mapping.
 
    [a,r,s,dim] = deal(argin{:});
		
    if any(isnan([r,s,dim]))        % optimize regularisation parameters
      defs = {0,0,[]};
      parmin_max = [1e-8,9.9999e-1;1e-8,9.9999e-1;1,size(a,2)];
      [w,r,s,dim] = regoptc(a,mfilename,{r,s,dim},defs,[3 2 1],parmin_max,testc([],'soft'),[1 1 0]);

    else % training

      islabtype(a,'crisp','soft'); % Assert A has the right labtype.
      isvaldfile(a,2,2); % at least 2 objects per class, 2 classes

      [m,k,c] = getsize(a);

      % If the subspace dimensionality is not given, set it to all dimensions.

      if (isempty(dim)), dim = k; end;

      dim = round(dim);
      if (dim < 1) | (dim > k)
        error ('Number of dimensions M should lie in the range [1,K].');
      end

      [U,G] = meancov(a);

      % Calculate means and priors.

      pars.mean  = +U;
      pars.prior = getprior(a);

      % Calculate class covariance matrices.

      pars.cov   = zeros(k,k,c);
      for j = 1:c
        F = G(:,:,j);

        % Regularize, if requested.

        if (s > 0) | (r > 0) 
          F = (1-r-s) * F + r * diag(diag(F)) +s*mean(diag(F))*eye(size(F,1));
        end

        % If DIM < K, extract the first DIM principal components and estimate
        % the noise outside the subspace.

        if (dim < k)
          dim = min(rank(F)-1,dim);	
          [eigvec,eigval] = preig(F); eigval = diag(eigval);
          [dummy,ind] = sort(-eigval);

          % Estimate sigma^2 as avg. eigenvalue outside subspace.
          sigma2 = mean(eigval(ind(dim+1:end)));

          % Subspace basis: first DIM eigenvectors * sqrt(eigenvalues).
          F = eigvec(:,ind(1:dim)) * diag(eigval(ind(1:dim))) * eigvec(:,ind(1:dim))' + ...
              sigma2 * eye(k);
        end

        pars.cov(:,:,j) = F;
      end

      w = normal_map(pars,getlab(U),k,c);
      w = setcost(w,a);

    end

    w = setname(w,mapname);
    
  end

return;
