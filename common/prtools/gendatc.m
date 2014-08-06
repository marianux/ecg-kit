%GENDATC Generation of two spherical classes with different variances
% 
% 	A = GENDATC(N,K,U,LABTYPE)
% 
% INPUT
%   N        Vector with class sizes (default: [50,50])
%   K        Dimensionality of the dataset (default: 2)
%   U        Mean of class 1 (default: 0)
%   LABTYPE  'crisp' or 'soft' labels (default: 'crisp')
%
% OUTPUT
%   A        Dataset
%
% DESCRIPTION
% Generation of a K-dimensional 2-class dataset A of N objects.
% Both classes are spherically Gaussian distributed.
%
% Class 1 has the identity matrix as covariance matrix and 
% mean U. If U is a scalar then [U,0,0,..] is used as class mean.
% Class 2 has also the identity matrix as covariance matrix, except
% for a variance of 4 for the first two features. Its mean is 0. 
% Class priors are P(1) = P(2) = 0.5.
%
% If N is a vector of sizes, exactly N(I) objects are generated
% for class I, I = 1,2.
%
% The default means result in a class overlap of 0.16.
% 
% LABTYPE defines the desired label type: 'crisp' or 'soft'. In the 
% latter case true posterior probabilities are set for the labels.
%
% Defaults: N = [50,50], K = 2, U = 0, LABTYPE = 'crisp'.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, PRDATASETS

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: gendatc.m,v 1.3 2009/01/27 13:01:42 duin Exp $

function A = gendatc(N,k,ma,labtype)

    if nargin < 1
    N = [50 50];
    prwarning(4,'size of classes not specified, assuming [50 50]');
  end
  if nargin < 2
    k=2;
    prwarning(4,'dimension not specified, assuming 2');
  end
  if nargin < 3
    ma=0;
    prwarning(4,'mean not specified, assuming 0');
  end
  if nargin < 4
    labtype = 'crisp';
    prwarning(4,'label type not specified, assuming crisp');
  end

  p = [0.5 0.5];
  N = genclass(N,p);
  % When a scalar ma is given, the mean vector should be generated:
  if (length(ma) == 1) & (k>1),
    ma=[ma,zeros(1,k-1)];
  end

  GA = eye(k);
  GB = eye(k); GB(1,1) = 9;
  if k > 1, GB(2,2) = 9; end
  mb = zeros(1,k);
  U = prdataset([ma;mb],[1 2]','prior',p);
  A = gendatgauss(N,U,cat(3,GA,GB));
  A = set(A,'name','Spherical Set');

  % Take care for the different types of labels:
  switch labtype
    case 'crisp'
         ;
    case 'soft'
      W = nbayesc(U,cat(3,GA,GB));
      targets = A*W*classc;
      A = setlabtype(A,'soft',targets);
    otherwise
      error(['Label type ' labtype ' not supported'])
    end
return
