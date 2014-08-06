%GENDATGAUSS (Formerly GAUSS) Generation of a multivariate Gaussian dataset
% 
% 	A = GENDATGAUSS(N,U,G,LABTYPE) 
%
% INPUT (in case of generation a 1-class dataset in K dimensions)
%   N		    Number of objects to be generated (default 50).
%   U		    Desired mean (vector of length K).
%   G       K x K covariance matrix. Default eye(K).
%   LABTYPE Label type (default 'crisp')
%
% INPUT (in case of generation a C-class dataset in K dimensions)
%   N       Vector of length C with numbers of objects per class.
%   U       C x K matrix with class means, or
%           Dataset with means, labels and priors of classes 
%           (default: zeros(C,K))
%   G       K x K x C covariance matrix of right size.
%           Default eye(K);
%   LABTYPE	Label type (default 'crisp')
%
% OUTPUT
%   A       Dataset containing multivariate Gaussian data
%
% DESCRIPTION
% Generation of N K-dimensional Gaussian distributed samples for C classes.
% The covariance matrices should be specified in G (size K*K*C) and the
% means, labels and prior probabilities can be defined by the dataset U with
% size (C*K). If U is not a dataset, it should be a C*K matrix and A will
% be a dataset with C classes.
%
% If N is a vector, exactly N(I) objects are generated for class I, 
% I = 1..C.
% 
% EXAMPLES
% 1. Generation of 100 points in 2D with mean [1 1] and default covariance
%    matrix: 
%
%        GENDATGAUSS(100,[0 0])
%
% 2. Generation of 50 points for each of two 1-dimensional distributions with
%    mean -1 and 1 and with variances 1 and 2:
%
%	     GENDATGAUSS([50 50],[-1;1],CAT(3,1,2))
%
%   Note that the two 1-dimensional class means should be given as a column
%   vector [1;-1], as [1 -1] defines a single 2-dimensional mean. Note that
%   the 1-dimensional covariance matrices degenerate to scalar variances,
%   but have still to be combined into a collection of square matrices using
%   the CAT(3,....) function.
%
% 3. Generation of 300 points for 3 classes with means [0 0], [0 1] and 
%    [1 1] and covariance matrices [2 1; 1 4], EYE(2) and EYE(2):
%
%      GENDATGAUSS(300,[0 0; 0 1; 1 1]*3,CAT(3,[2 1; 1 4],EYE(2),EYE(2)))
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>) 
% DATASETS, PRDATASETS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function a = gendatgauss(n,u,g,labtype)

		if (nargin < 1)
		prwarning (2,'number of samples not specified, assuming N = 50'); 	
		n = 50;
		end
		cn = length(n);
	if (nargin < 2)
		prwarning (2,'means not specified; assuming one dimension, mean zero');
		u = zeros(cn,1); 
	end;
	if (nargin < 3)
		prwarning (2,'covariances not specified, assuming unity');
	 	g = eye(size(u,2)); 
	end
	if (nargin < 4)
		prwarning (3,'label type not specified, assuming crisp');
		labtype = 'crisp'; 
	end

	% Return an empty dataset if the number of samples requested is 0.

	if (length(n) == 1) & (n == 0)
		a = prdataset([]); 
		return
	end

	% Find C, desired number of classes based on U and K, the number of 
	% dimensions. Make sure U is a dataset containing the means.

	if (isa(u,'prdataset'))
		[m,k,c] = getsize(u);
		lablist = getlablist(u);
		p = getprior(u);
		if c == 0
			u = double(u);
		end
	end
	if isa(u,'double')
		[m,k] = size(u); 		
		c = m;
		lablist = genlab(ones(c,1));
		u = prdataset(u,lablist);
		p = ones(1,c)/c;
	end

	if (cn ~= c) & (cn ~= 1)
		error('The number of classes specified by N and U does not match');
	end

 	% Generate a class frequency distribution according to the desired priors.
	n = genclass(n,p);

	% Find CG, the number of classes according to G. 
	% Make sure G is not a dataset.

	if (isempty(g))
		g = eye(k); 
		cg = 1;
	else
		g = real(+g); [k1,k2,cg] = size(g);
		if (k1 ~= k) | (k2 ~= k)
			error('The number of dimensions of the means U and covariance matrices G do not match');
		end
		if (cg ~= m & cg ~= 1)
			error('The number of classes specified by the means U and covariance matrices G do not match');
		end
	end

	% Create the data A by rotating and scaling standard normal distributions 
	% using the eigenvectors of the specified covariance matrices, and adding
	% the means.

	a = [];
	for i = 1:m
		j = min(i,cg);						% Just in case CG = 1 (if G was not specified).

		% Sanity check: user can pass non-positive definite G.		
		[V,D] = preig(g(:,:,j)); V = real(V); D = real(D); D = max(D,0);
		a = [a; randn(n(i),k)*sqrt(D)*V' + repmat(+u(i,:),n(i),1)];
	end

	% Convert A to dataset by adding labels and priors.

	labels = genlab(n,lablist);
	a = prdataset(a,labels,'lablist',lablist,'prior',p);

	% If non-crisp labels are requested, use output of Bayes classifier.
	switch (labtype)
		case 'crisp'
			;
		case 'soft'
			w = nbayesc(u,g); 		
			targets = a*w*classc;
			a = setlabtype(a,'soft',targets);
		otherwise
			error(['Label type ' labtype ' not supported'])
	end

	a = setname(a,'Gaussian Data');

return
