%TESTN Error estimate of discriminant for normal distribution.
% 
% 	E = TESTN(W,U,G,N)
% 
% INPUT
%   W  Trained classifier mapping
%   U  C x K dataset with C class means, labels and priors (default: [0 .. 0])
%   G  K x K x C matrix with C class covariance matrices (default: identity)
%   N  Number of test examples (default 10000)
%
% OUTPUT
%   E  Estimated error
%      
% DESCRIPTION
% This routine estimates as good as possible the classification error
% of Gaussian distributed problems with known means and covariances. N
% normally distributed data vectors with means, labels and prior
% probabilities defined by the dataset U (size [C,K]) and covariance
% matrices G (size [K,K,C]) are generated with the specified labels and are
% tested against the discriminant W. The fraction of incorrectly classified
% data vectors is returned. If W is a linear 2-class discriminant and N 
% is not specified, the error is computed analytically.
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, QDC, NBAYESC, TESTC

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: testn.m,v 1.3 2009/01/27 13:01:42 duin Exp $

function e = testn(w,U,G,m)

		
	% Assert that W is a trained mapping.
	istrained(w);	[k,c] = size(w);

	if (nargin < 4) & ~( isaffine(w) & (c == 2) )
		prwarning(4,'number of samples N to test with not specified, assuming 10000');
		m = 10000; 
	end
	if (nargin < 3)
		prwarning(3,'covariance matrices not specified, assuming identity for each class');
		G = eye(k); 
	end
	if (nargin < 2)
		error('Class means are not specified');
	else
		[cc,k] = size(U); p = getprior(U); u = +U;
		if (cc ~= c)
			error('Specified number of class means does not fit classifier')
		end
	end

	% If a single covariance is specified, use it for each class.

	if (length(size(G)) == 2)
		g = G;
		for j = 2:c
			G = cat(3,G,g);
		end
	else
		if (size(G,3) ~= c)
			error('Specified number of covariances does not fit classifier')
		end
	end
	
	% Check for analytical case: linear classifier and no data samples requested.

	if (nargin < 4) & (isaffine(w)) & (c == 2)

		prwarning(2,'computing Bayes error analytically');

		e = 0;
		v = w.data.rot(:,1); f = w.data.offset(1);
		for j = 1:c

			% Find the index J of mean j stored in U.
			J = findnlab(U,j);
			if (length(J) ~= 1)
				error('The mean U does not contain correct labels.')
			end

			q = real(sqrt(v'*G(:,:,j)*v)); % std dev in direction of classifier
			d = (2*j-3)*(v'*u(J,:)'+f);    % normalized distance to classifier
			if (q == 0)
				if (d >= 0), e = e + p(j); end
			else
				e = e + p(j) * (erf(d/(q*sqrt(2)))/2 + 0.5);
			end

		end
	else 

		% Cannot solve analytically: generate data and test.
		a = gendatgauss(m,U,G); 		% Generate normally distributed data.
		a = setlablist(a,getlab(w));
		e = testc(a,w);   		% Find the error of discriminant on the data.

	end

return
