%GENTRUNK Generation of Trunk's classification problem of 2 Gaussian classes
% 
%   A = GENTRUNK(N,K)
% 
% INPUT
%   N       Dataset size, or 2-element array of class sizes (default: [50 50]).
%   K       Dimensionality of the dataset to be generated (default: 10).
%
% OUTPUT
%   A       Dataset.
%
% DESCRIPTION
% Generation of a K-dimensional 2-class dataset A of N objects. Both classes 
% are Gaussian distributed with the idenity matrix as covariance matrix.
% The means of the first class are defined by ua(j) = 1/sqrt(j). The means
% for the second class are ub = -ua. These means are such that the Nearest
% Mean Classifier always shows peaking for a finite training set.
%
% REFERENCES
% 1. G.V. Trunk, A Problem of Dimensionality: A Simple Example, IEEE Trans. 
% Pattern Analysis and Machine Intelligence, vol. 1, pp. 306-307, 1979
% 2. A.K. Jain, R.P.W. Duin, and J. Mao, Statistical Pattern Recognition: 
% A Review, IEEE Transactions on Pattern Analysis and Machine Intelligence, 
% vol. 22, pp. 4-37, 2000.
%
% EXAMPLE
% a = gentrunk([1000 1000],200);
% e = clevalf(a,nmc,[1:9 10:5:25 50:25:200],[5 5],25);
% plote(e)
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, PRDATASETS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: gentrunk.m,v 1.2 2009/01/27 13:01:42 duin Exp $

function A = gendats (N,k)

		if (nargin < 1), N = [50 50]; end
	if (nargin < 2), k = 10;	end

	% Set equal priors and generate random class sizes according to these.
	p = [0.5 0.5]; N = genclass(N,p);	

	% Unit covariance matrices
	% GA = eye(k); GB = eye(k);
  % Trunk means
  ma = 1./sqrt(1:k);
  mb = -ma;
  
  A = prdataset([randn(N(1),k) + repmat(ma,N(1),1); randn(N(2),k) + repmat(mb,N(2),1)]);
  A = prdataset(A,genlab(N));
  A = setprior(A,p);
  
% 	U = prdataset([ma;mb],[1 2]');
% 	U = setprior(U,p);
% 
% 	% Create dataset.
% 	A = gendatgauss(N,U,cat(3,GA,GB));
	A = setname(A,'Trunk''s Problem');

return
