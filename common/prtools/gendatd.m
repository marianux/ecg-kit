%GENDATD Generation of 'difficult' normally distributed classes
% 
%   A = GENDATD(N,K,D1,D2,LABTYPE)
%
% INPUT
%   N        Number of objects in each of the classes (default: [50 50])
%   K        Dimensionality of the dataset (default: 2)
%   D1       Difference in mean in feature 1 (default: 3)
%   D2       Difference in mean in feature 2 (default: 3)
%   LABTYPE  'crisp' or 'soft' labels (default: 'crisp').
%
% OUTPUT
%   A        Generated dataset
%
% DESCRIPTION
% Generation of a K-dimensional 2-class dataset A of N objects.
% Class variances are very different for the first two dimensions.
% Separation is thereby, for small sample sizes, 'difficult'. 
% 
% D1 is the difference between the means for the first feature, D2
% is the difference between the means for the second feature. In all
% other directions the means are equal. The two covariance matrices
% are equal with a variance of 1 in all directions except for the
% second feature, which has a variance of 40. The first two feature
% are rotated over 45 degrees to construct a strong correlation.
% Class priors are P(1) = P(2) = 0.5.
%
% If N is a vector of sizes, exactly N(I) objects are generated
% for class I, I = 1,2.
%
% LABTYPE defines the desired label type: 'crisp' or 'soft'. In the 
% latter case true posterior probabilities are set for the labels.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, PRDATASETS

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: gendatd.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function A = gendatd(N,k,d1,d2,labtype)

		
	if nargin < 5, labtype = 'crisp'; end
	if nargin < 4, d2 = 3; end
	if nargin < 3, d1 = 3; end
	if nargin < 2,  k = 2; end
	if nargin < 1, N = [50 50]; end

	if k < 2,
		error('Number of features should be larger than 1'),
	end
	V = ones(1,k); V(2) = 40; V = sqrt(V);
	p = [0.5 0.5];
	N = genclass(N,p);	
	ma = zeros(1,k);
	mb = zeros(1,k); mb(1:2) = [d1, d2];
	A = [randn(N(1),k).*V(ones(1,N(1)),:) + ma(ones(1,N(1)),:); ...
	     randn(N(2),k).*V(ones(1,N(2)),:) + mb(ones(1,N(2)),:)];
	A(:,1:2) = A(:,1:2)*[1 -1; 1 1]./sqrt(2);
	lab = genlab(N);
	A = prdataset(A,lab,'name','Difficult Dataset','prior',p);

	switch labtype
	 case 'crisp'
	  ;
	 case 'soft'
	  U = prdataset([ma(1:2);mb(1:2)],getlablist(A));
	  G = diag(V(1:2));
	  W = nbayesc(U,G);
	  targets = A(:,1:2)*W*classc;
	  A = setlabtype(A,'soft',targets);
	 otherwise
	  error(['Label type ' labtype ' not supported'])
	end

	return
