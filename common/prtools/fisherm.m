%FISHERM Optimal discrimination linear mapping (Fisher mapping, LDA)
%
%  W = FISHERM(A,N,ALF)
% 
% INPUT
%  A   Dataset
%  N   Number of dimensions to map to, N < C, where C is the number of classes
%       (default: min(C,K)-1, where K is the number of features in A)
%  ALF Preserved variance in the pre-whitening step
%
% OUTPUT
%  W   Fisher mapping
%
% DESCRIPTION  
% Finds a mapping of the labeled dataset A onto an N-dimensional linear
% subspace such that it maximizes the the between scatter over the within
% scatter (also called the Fisher mapping [1 or LDA]). Note that N should be 
% less than the number of classes in A. If supplied, ALF determines the 
% preserved variance in the prewhitening step (i.e. removal of insignificant 
% eigenvectors in the within-scatter, the EFLD procedure [2]), see KLMS.
%
% The resulting mapping is not orthogonal. It may be orthogonalised by ORTH.
% 
% REFERENCES
% [1] K. Fukunaga, Introduction to statistical pattern recognition, 2nd
% ed., Academic Press, New York, 1990.
% [2] C. Liu and H. Wechsler, Robust Coding Schemes for Indexing and Retrieval
% from Large Face Databases, IEEE Transactions on Image Processing, vol. 9, 
% no. 1, 2000, 132-136.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, NLFISHERM, KLM, PCAM, KLMS, ORTH

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: fisherm.m,v 1.4 2007/04/21 23:05:59 duin Exp $

function w = fisherm(a,n,alf)

	    
	if (nargin < 3)
		alf = [];
		prwarning (4, 'no pre-whitening reduction specified, so take all features');
	end
	if (nargin < 2)
		prwarning (4, 'number of dimensions to map to not specified, assuming min(C,K)-1');
		n = []; 
	end

	% If no arguments are specified, return an untrained mapping.

	if (nargin < 1) | (isempty(a))
		w = prmapping('fisherm',{n,alf});
		w = setname(w,'Fisher mapping');
		return
	end

	islabtype(a,'crisp','soft');
	a = testdatasize(a);
	isvaldset(a,2,2); % at least 2 objects per class, 2 classes

	[m,k,c] = getsize(a);
	% If N is not given, set it. 
	%DXD   Question: why do we subtract here 1 from k? Should it not be
	%      the next one??
	%if (isempty(n)), n = min(k,c)-1; end
	if (isempty(n)), n = min(k,c-1); end
	if (n >= m) | (n >= c) | (n > k)
		error('The dataset is too small or has too few classes for the requested dimensionality')
	end

	% A Fisher mapping is determined by the eigenvectors of inv(W)*B, where W
  % is the within-scatter matrix, i.e. the average covariance matrix
  % weighted by the prior probabilities, and B is the between-scatter
  % matrix. 

	% To simplify computations, W can be set to the identity matrix if A is 
	% centered and sphered by KLMS.
	
	a = setprior(a,getprior(a)); % set priors to avoid unnecessary warnings
	v = klms(a,alf); a = a*v;			
	
	% Calculate eigenvectors of B.

	u = pcam(meancov(a),n);  
	
	if (n == 0)								
		%DXD: Careful here: I wanted to make the option n==0 behave the
	   %     same as in pcam(x,0) and bhatm(x,0). Unfortunately, it was
		%     already defined. So I overrule it now!!!! Please correct
		%     if you don't agree!
		w = u;
		return
		%w = v;						% Only one feature is available	(k=1)
	else
		w = v*u;					% Construct final mapping.
	end
	w = setname(w,'Fisher mapping');

return
