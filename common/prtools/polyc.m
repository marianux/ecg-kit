%POLYC Polynomial Classification
% 
%   W = polyc(A,CLASSF,N,S)
% 
% INPUT
%   A       Dataset
%   CLASSF  Untrained classifier (optional; default: FISHERC)
%   N       Degree of polynomial (optional; default: 1)
%   S       1/0, 1 indicates that 2nd order combination terms should be used
%           (optional; default: 0)
%
% OUTPUT
%	  W       Trained classifier
%
% DESCRIPTION
% Adds polynomial features to the dataset A and runs the untrained 
% classifier CLASSF. Combinations of 2nd order terms may be constructed.
% For higher order terms no combinations are generated.
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: polyc.m,v 1.3 2007/04/13 09:30:54 duin Exp $

function W = polyc(a,classf,n,s)

		if nargin < 4, 
		prwarning(4,'S is not specified, assuming 0, i.e. no higher order terms.'); 	
		s = 0; 
	end
	if nargin < 3 | isempty(n), 
		prwarning(4,'N, the degree of polynomial is not specified, assuming 1.'); 	
		n = 1; 
	end
	if (nargin < 2) | (isempty(classf)), 
		prwarning(4,'CLASSF is not specified, assuming FISHERC.'); 	
		classf = fisherc; 
	end

	% No data, return an untrained mapping.
	if nargin < 1 | isempty(a)
		W = prmapping('polyc','untrained',{classf,n,s});
		W = setname(W,'Polynomial classifier');
		return;
	end

	% Checks for the classifier.

	ismapping(classf);
	isuntrained(classf);
	% Construct the 'feature-construction'-matrix such that:
	% 		x_new = x*cmapm(k,P);
	% See CMAPM.

	isvaldfile(a,1,2); % at least 1 object per class, 2 classes
	a = testdatasize(a,'features');
	a = setprior(a,getprior(a));

	[m,k] = size(a);
	
	% First all the terms involving the same degree are defined.

	% Degree 1
	P = eye(k);
	% Degree 2 and higher.
	for j = 2:n
		P = [P; j*eye(k)];
	end

	% On request, all the terms involving combinations of different degrees.
	if s & (k > 1) & n > 1
		% All pairs of features (including their higher degree variants) have
		% to be included.
		Q = zeros(k*(k-1)/2,k);
		n = 0;

		for j1 = 1:k
			for j2 = j1+1:k
				n = n+1;
				Q(n,j2) = Q(n,j2)+1;
				Q(n,j1) = Q(n,j1)+1;
			end
		end
		P = [P;Q];
	end
	
	% Define the mapping from P, apply it to A and train the classifier.
	v = cmapm(k,P);
	w = a*v*classf;
	W = v*w;
	W = setname(W,'Polynomial classifier');
	W = setcost(W,a);

return;
