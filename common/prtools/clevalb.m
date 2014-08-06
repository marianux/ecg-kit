%CLEVALB Classifier evaluation (learning curve), bootstrap version
% 
% 	E = CLEVALB(A,CLASSF,TRAINSIZES,NREPS)
%
% INPUT
%   A           Training dataset
%   CLASSF      Classifier to evaluate
%   TRAINSIZES  Vector of class sizes, used to generate subsets of A
%               (default [2,3,5,7,10,15,20,30,50,70,100])
%   NREPS       Number of repetitions (default 1)
%
% OUTPUT
%   E           Error structure (see PLOTE)
%
% DESCRIPTION 
% Generates at random, for all class sizes defined in TRAINSIZES, training
% sets out of the dataset A and uses these for training the untrained
% classifier CLASSF. CLASSF may also be a cell array of untrained
% classifiers; in this case the routine will be run for all of them. The
% resulting trained classifiers are tested on all objects in A. This
% procedure is then repeated N times.
%
% Training set generation is done "with replacement" and such that for each
% run the larger training sets include the smaller ones and that for all
% classifiers the same training sets are used.
% 
% If CLASSF is fully deterministic, this function uses the RAND random
% generator and thereby reproduces if its seed is reset (see RAND). 
% If CLASSF uses RANDN, its seed may have to be set as well.
%
% Use FID = 1 to report progress to the command window.
% 
% EXAMPLES
% See PREX_CLEVAL.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, CLEVALB, TESTC, PLOTE

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: clevalb.m,v 1.4 2008/07/03 09:05:50 duin Exp $

function e = clevalb(a,classf,learnsizes,nreps,fid) 
  
  % use of fid is outdated
  if (nargin < 4) || isempty(nreps)
		prwarning(2,'number of repetitions not specified, assuming NREPS = 1');
		nreps = 1; 
	end;
  if (nargin < 3) || isempty(learnsizes)
		prwarning(2,'vector of training set class sizes not specified, assuming [2,3,5,7,10,15,20,30,50,70,100]');
 		learnsizes = [2,3,5,7,10,15,20,30,50,70,100]; 
	end;

	% If a single mapping is given, convert it to a 1 x 1 cell array.

  if (ismapping(classf)), classf = {classf}; end

	% Correct for old argument order.

  if (isdataset(classf)) & (ismapping(a)) 
  	tmp = a; a = classf; classf = {tmp};
  end
  if (isdataset(classf)) & (iscell(a)) & (ismapping(a{1})) 
  	tmp = a; a = classf; classf = tmp;
  end
  if ~iscell(classf), classf = {classf}; end

	% Assert that all is right.
  isdataset(a); ismapping(classf{1});

	% Remove requested class sizes that are larger than the size of the
	% smallest class.

  mc = classsizes(a); [m,k,c] = getsize(a);
	toolarge = find(learnsizes >= min(mc));
	if (~isempty(toolarge))
		prwarning(2,['training set class sizes ' num2str(learnsizes(toolarge)) ...
								 ' larger than the minimal class size in A; removed them']);
	  learnsizes(toolarge) = [];
	end
  learnsizes = learnsizes(:)';

	% Fill the error structure.

  nw = length(classf(:));
  datname = getname(a);

  e.n       = nreps;
  e.error   = zeros(nw,length(learnsizes));
  e.std     = zeros(nw,length(learnsizes));
  e.xvalues = learnsizes(:)';
  e.names   = [];
  e.xlabel  = 'Training set size';
  if (nreps > 1)
  	e.ylabel= ['Averaged error (' num2str(nreps) ' experiments)'];
  elseif (nreps == 1)
  	e.ylabel = 'Error';
  else
		error('Number of repetitions NREPS should be >= 1.');
	end;
  if (~isempty(datname))
  	e.title = ['Bootstrapped learning curve on ' datname];
  end
  if (learnsizes(end)/learnsizes(1) > 20)
  	e.plot = 'semilogx'; 				% If range too large, use a log-plot for X.
  end
  
	% Report progress.
	
	s1 = sprintf('clevalb: %i classifiers: ',nw);
	prwaitbar(nw,s1);

	% Store the seed, to reset the random generator later for different
	% classifiers.

  seed = randreset;

	% Loop over all classifiers (with index WI).

  for wi = 1:nw

  	isuntrained(classf{wi}); 
    name = getname(classf{wi});
    e.names = char(e.names,name);
    prwaitbar(nw,wi,[s1 name]);

		% E1 will contain the error estimates.

  	e1 = zeros(nreps,length(learnsizes));

		% Take care that classifiers use same training set.

  	randreset(seed); seed2 = seed;

		% For NREPS repetitions...
		s2 = sprintf('cleval: %i repetitions: ',nreps);
		prwaitbar(nreps,s2);

  	for i = 1:nreps

			prwaitbar(nreps,i,[s2 int2str(i)]);
			% Store the randomly permuted indices of samples of class CI to use in
			% this training set in JR(CI,:).

  		JR = zeros(c,max(learnsizes));
  		for ci = 1:c
				JC = findnlab(a,ci);
				
				% Necessary for reproducable training sets: set the seed and store
				% it after generation, so that next time we will use the previous one.
  			randreset(seed2);		

  			R = ceil(rand(1,max(learnsizes))*length(JC));
  			JR(ci,:) = JC(R)';

  			seed2 = randreset; 
  		end

  		li = 0;										% Index of training set.
			nlearns = length(learnsizes);
			s3 = sprintf('cleval: %i sizes: ',nlearns);
			prwaitbar(nreps,s3);
  		for j = 1:nlearns
				nj = learnsizes(j);
        prwaitbar(nlearns,j,[s3 int2str(j) ' (' int2str(nj) ')']);
        li = li + 1; 

				% J will contain the indices for this training set.

				J = [];				
  			for ci = 1:c
  				J = [J;JR(ci,1:nj)']; 
  			end;
				
				% Train classifier CLASSF{WI} on this training set and calculate
				% error.
				
  			W = a(J,:)*classf{wi}; 
    		e1(i,li) = testc(a,W);

  		end
      prwaitbar(0);

  	end
    prwaitbar(0);

		% Calculate average error and standard deviation for this classifier
		% (or set the latter to zero if there's been just 1 repetition).

  	e.error(wi,:) = mean(e1,1);
  	if (nreps == 1)
  		e.std(wi,:) = zeros(1,size(e.std,2));
  	else
  		e.std(wi,:) = std(e1)/sqrt(nreps);
  	end

  end
  prwaitbar(0);
  
  % The first element is the empty string [], remove it.
  e.names(1,:) = [];

return
