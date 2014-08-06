%CLEVALF Classifier evaluation (feature size curve)
% 
%   E = CLEVALF(A,CLASSF,FEATSIZES,TRAINSIZE,NREPS,S,TESTFUN)
% 
% INPUT
%   A          Training dataset.
%   CLASSF     The untrained classifier to be tested.
%   FEATSIZES  Vector of feature sizes (default: all sizes)
%   TRAINSIZE  Number of objects/fraction of training set size (see GENDAT)
%              or generator mapping.
%   NREPS      Number of repetitions (default: 1)
%   S          Independent test dataset (optional)
%   TESTFUN    Mapping,evaluation function (default classification error)
%
% OUTPUT
%   E          Structure with results
%              See PLOTE for a description
%
% DESCRIPTION
% Generates at random for all feature sizes stored in FEATSIZES training
% sets of the given TRAINSIZE out of the dataset A. See GENDAT for the
% interpretation of TRAINSIZE. These are used for training the untrained
% classifier CLASSF. The result is tested by all unused ojects of A, or,
% if given, by the test dataset S. This is repeated N times. If no testset
% is given and if LEARNSIZE is not given or empty, the training set is
% bootstrapped. If a testset is given, the default training set size is 
% the entire training set. Default FEATSIZES: all feature sizes. 
% The mean erors are stored in E. The observed standard deviations are 
% stored in S. The default test routine is classification error estimation 
% by TESTC([],'crisp'). 
%
% See CLEVALFS for how to construct feature curves in addition with
% automatic feature extraction / selection.
% 
% This function uses the RAND random generator and thereby reproduces only
% if its seed is saved and reset.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>) 
% MAPPINGS, DATASETS, CLEVAL, CLEVALFS, TESTC, PLOTE, GENDAT

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: clevalf.m,v 1.8 2008/07/03 09:05:50 duin Exp $

function e = clevalf(a,classf,featsizes,learnsize,n,Tset,testfun)
	
		[m,k] = size(a);
  if (nargin < 7) | isempty(testfun)
    testfun = testc([],'crisp');
  end;
	if nargin < 6, Tset = []; end
	if nargin < 5, n = 1; end;
	if nargin < 4 | isempty(learnsize), learnsize = [0.5]; end
	if nargin < 3 | isempty(featsizes), featsizes = [1:k]; end
	
	if isdataset(classf) & ismapping(a) % correct for old order
		dd = a; a = classf; classf = {dd};
	end
	if isdataset(classf) & iscell(a) & ismapping(a{1}) % correct for old order
		dd = a; a = classf; classf = dd;
	end
	
	if ~iscell(classf), classf = {classf}; end
	isdataset(a);
	ismapping(classf{1});
  a = setprior(a,getprior(a));

	if ~isempty(Tset), isdataset(Tset); T = Tset; end

	[m,k,c] = getsize(a);
	featsizes(find(featsizes > k)) = [];
	featsizes = featsizes(:)';

	if ~ismapping(learnsize) && (length(learnsize) > 1 & length(learnsize) ~= c)
		error('Learnsize should be scalar or a vector with length equal to the class size')
	end
	
	r = length(classf(:));
	e.error = zeros(r,length(featsizes));
	e.std   = zeros(r,length(featsizes));
	e.xvalues = featsizes;
	e.n = n;
	datname = getname(a);
	if ~isempty(datname)
		e.title = ['Feature curve for ' getname(a)];
	end
	e.xlabel= 'Feature size';
	if n > 1
		e.ylabel= ['Averaged error (' num2str(n) ' experiments)'];
	else
		e.ylabel = 'Error';
	end

	if featsizes(end)/featsizes(1) > 20
		e.plot = 'semilogx';
  else
    e.plot = 'plot';
	end
	e.names = [];

	s1 = sprintf('clevalf: %i classifiers: ',r);
	prwaitbar(r,s1);
	
	e1 = zeros(n,length(featsizes));
	seed = randreset;

	% loop over all classifiers
	
	for q = 1:r
		isuntrained(classf{q});
		name = getname(classf{q});
		prwaitbar(r,q,[s1 name]);
		e.names = char(e.names,name);
		e1 = zeros(n,length(featsizes));
		randreset(seed);  % take care that classifiers use same training set
		seed2 = seed;
		s2 = sprintf('clevalf: %i repetitions: ',n);
		prwaitbar(n,s2);
		for i = 1:n
			prwaitbar(n,i,[s2 int2str(i)]);
			randreset(seed2);
			if isempty(Tset)
				[b,T] = gendat(a,learnsize);
			elseif ~isempty(learnsize)
				b = gendat(a,learnsize);
			else
				b = a;
			end
			seed2 = randreset;
			nfeatsizes = length(featsizes);
			s3 = sprintf('clevalf: %i feature sizes: ',nfeatsizes);
			prwaitbar(nfeatsizes,s2);
			for j=1:nfeatsizes
				f = featsizes(j);
				prwaitbar(nfeatsizes,j,[s3 int2str(j) ' (' int2str(f) ')']);
				e1(i,j) = T(:,1:f)*(b(:,1:f)*classf{q})*testfun;
			end
			prwaitbar(0)
		end
		prwaitbar(0)
		e.error(q,:) = mean(e1,1);
		if n == 1
			e.std(q,:) = zeros(1,size(e.std,2));
		else
			e.std(q,:) = std(e1)/sqrt(n);
		end
	end
	prwaitbar(0)

	e.names(1,:) = [];
	
	return
	
