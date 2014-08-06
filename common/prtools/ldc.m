%LDC Linear Bayes Normal Classifier (BayesNormal_1)
%
%  [W,R,S,M] = LDC(A,R,S,M)
%  [W,R,S,M] = A*LDC([],R,S,M);
%  [W,R,S,M] = A*LDC(R,S,M);
% 
% INPUT
%  A    Dataset
%  R,S  Regularization parameters, 0 <= R,S <= 1
%      (optional; default: no regularization, i.e. R,S = 0)
%  M    Dimension of subspace structure in covariance matrix (default: K,
%       all dimensions)
%
% OUTPUT
%  W  Linear Bayes Normal Classifier mapping
%  R  Value of regularization parameter R as used 
%  S  Value of regularization parameter S as used
%  M  Value of regularization parameter M as usedd
%
% DESCRIPTION	
% Computation of the linear classifier between the classes of the dataset A
% by assuming normal densities with equal covariance matrices. The joint
% covariance matrix is the weighted (by a priori probabilities) average of
% the class covariance matrices. R and S (0 <= R,S <= 1) are regularization
% parameters used for finding the covariance matrix G by:
%
%  G = (1-R-S)*G + R*diag(diag(G)) + S*mean(diag(G))*eye(size(G,1))
%
% This covariance matrix is then decomposed as 
%
%  G = W*W' + sigma^2 * eye(K)
%
% where W is a K x M matrix containing the M leading principal components
% and sigma^2 is the mean of the K-M smallest eigenvalues. The use of soft labels 
% is supported. The classification A*W is computed by NORMAL_MAP.
%
% If R, S or M is NaN the regularisation parameter is optimised by REGOPTC.
% The best result are usually obtained by R = 0, S = NaN, M = [], or by
% R = 0, S = 0, M = NaN (which is for problems of moderate or low dimensionality
% faster). If no regularisation is supplied a pseudo-inverse of the
% covariance matrix is used in case it is close to singular.
%
% Note that A*(KLMS([],N)*NMC) performs a similar operation by first
% pre-whitening the data in an N-dimensional space, followed by the
% nearest mean classifier. The regularization controlled by N is different
% from the above in LDC as it entirely removes small variance directions.
%
% To some extend LDC is also similar to FISHERC.
%
% EXAMPLES
% See PREX_PLOTC.
% a = gendatd;  % generate Gaussian distributed data in two classes
% w = ldc(a);   % compute a linear classifier between the classes
% scatterd(a);  % make a scatterplot
% plotc(w)      % plot the classifier
%
% REFERENCES
% 1. R.O. Duda, P.E. Hart, and D.G. Stork, Pattern classification, 2nd edition, 
%    John Wiley and Sons, New York, 2001.
% 2. A. Webb, Statistical Pattern Recognition, John Wiley & Sons, New York, 2002.
% 3. C. Liu and H. Wechsler, Robust Coding Schemes for Indexing and Retrieval
%    from Large Face Databases, IEEE Transactions on Image Processing, vol. 9, 
%    no. 1, 2000, 132-136.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, REGOPTC, NMC, NMSC, LDC, UDC, QUADRC, NORMAL_MAP, FISHERC

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: ldc.m,v 1.11 2010/02/08 15:31:48 duin Exp $

function [W,r,s,dim] = ldc(varargin)

  mapname = 'Bayes-Normal-1';
  argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],0,0,[]);
  
  if mapping_task(argin,'definition')
    
    W = define_mapping(argin,'untrained',mapname);
    
	elseif mapping_task(argin,'training')			% Train a mapping.
 
    [a,r,s,dim] = deal(argin{:});
	
    if any(isnan([r,s,dim]))		% optimize regularisation parameters
      defs = {0,0,[]};
      parmin_max = [1e-8,9.9999e-1;1e-8,9.9999e-1;1,size(a,2)];
      [W,r,s,dim] = regoptc(a,mfilename,{r,s,dim},defs,[3 2 1],parmin_max,testc([],'soft'),[1 1 0]);

    else % training with known reg pars

      islabtype(a,'crisp','soft');
      isvaldfile(a,2,2); % at least 2 object per class, 2 classes

      [m,k,c] = getsize(a);

      % Calculate mean vectors, priors and the covariance matrix G.
      %DXD: for large datasets (high dim, many classes) it is impossible
      %to keep all cov. matrices into memory. We have to loop over the
      %individual classes:
      %[U,G] = meancov(a);
      U = meancov(a);
      w.mean	= +U;
      w.prior = getprior(a);
      %G = reshape(sum(reshape(G,k*k,c)*w.prior',2),k,k);
      [tmpU,G] = meancov(seldat(a,1)); 
      G = w.prior(1)*G;
      for i=2:c
        [tmpU,tmpG] = meancov(seldat(a,i)); 
        G = G + w.prior(i)*tmpG;
      end
      clear tmpG;

      % Regularize 
      if (s > 0) | (r > 0)
        G = (1-r-s)*G + r * diag(diag(G)) + s*mean(diag(G))*eye(size(G,1));
      end	

      if (dim < k)
        dim = min(rank(G)-1,dim);
        [eigvec,eigval] = preig(G); eigval = diag(eigval);
        [dummy,ind] = sort(-eigval);
        % Estimate sigma^2 as avg. eigenvalue outside subspace.
        sigma2 = mean(eigval(ind(dim+1:end)));
        % Subspace basis: first DIM eigenvectors * sqrt(eigenvalues).
        G = eigvec(:,ind(1:dim)) * diag(eigval(ind(1:dim))) * eigvec(:,ind(1:dim))' + ...
            sigma2 * eye(k);
      end

      w.cov = G;
      W = normal_map(w,getlab(U),k,c);
      W = setcost(W,a);

    end

    W = setname(W,mapname);
    
  end
	
return
