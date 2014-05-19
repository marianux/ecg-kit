%CLEVALS Classifier evaluation (feature size/learning curve), bootstrap possible
% 
% 	E = CLEVALS(A,CLASSF,FEATSIZE,TRAINSIZES,NREPS,T)
%
% INPUT
%   A          Training dataset
%   CLASSF     Classifier to evaluate
%   FEATSIZE   Vector of feature sizes
%                (default: 1:K, where K is the number of features in A)
%   TRAINSIZES Vector of class sizes, used to generate subsets of A
% 	           	 (default [2,3,5,7,10,15,20,30,50,70,100])
%   NREPS 		 Number of repetitions (default 1)
%   T      		 Tuning set, or 'bootstrap' (default [], i.e. use remaining
%             		objects in A)
%
% OUTPUT
%   E       	 Error structure (see PLOTE)
%
% DESCRIPTION 
% Generates at random, for all feature sizes defined in FEATSIZES or all
% class sizes defined in TRAINSIZES, training sets out of the dataset A and
% uses these for training the untrained classifier CLASSF. CLASSF may also
% be a cell array of untrained classifiers; in this case the routine will be
% run for all of them. The resulting trained classifiers are tested on all
% objects in A. This procedure is then repeated N times.
%
% Training set generation is done "with replacement" and such that for each
% run the larger training sets include the smaller ones and that for all
% classifiers the same training sets are used.
% 
% If CLASSF is fully deterministic, this function uses the RAND random
% generator and thereby reproduces if its seed is reset (see RAND). 
% If CLASSF uses RANDN, its seed may have to be set as well.
% 
% EXAMPLES
% See PREX_CLEVAL.
%
% SEE ALSO
% MAPPINGS, DATASETS, CLEVALB, TESTC, PLOTE

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: clevals.m,v 1.4 2008/07/03 09:05:50 duin Exp $

function e = clevals(a,classf,featsizes,learnsizes,nreps,t,fid) 

	  
  prwaitbar

  % use of fid is outdated
	if (nargin < 6)
		prwarning(2,'no tuning set supplied, bootstrapping');
		t = [];
	end;
  if (nargin < 5)
		prwarning(2,'number of repetitions not specified, assuming NREPS = 1');
		nreps = 1; 
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

	% Defaults for size arrays.

	if (nargin < 4) 
		prwarning(2,'vector of training set class sizes not specified, assuming [2,3,5,7,10,15,20,30,50,70,100]');
 		learnsizes = [2,3,5,7,10,15,20,30,50,70,100]; 
	end;
  if (nargin < 3)
		prwarning(2,'vector of feature sizes not specified, assuming K');
 		featsizes = k;
	end;

	learncurve = 0; featcurve = 0;
	if (isempty(featsizes))									% Learning curve.
		toolarge = find(learnsizes >= min(mc));
		if (~isempty(toolarge))
			prwarning(2,['training set class sizes ' num2str(learnsizes(toolarge)) ...
									 ' larger than the minimal class size in A; removed them']);
		  learnsizes(toolarge) = [];
		end
		learnsizes = learnsizes(:)';
		featsizes  = k;
		sizes = learnsizes;
		learncurve = 1;
	else																		% Feature size curve.
		toolarge = find(featsizes > k);
		if (~isempty(toolarge))
			prwarning(2,['feature sizes ' num2str(featsizes(toolarge)) ...
									 ' larger than number of features in A; removed them']);
		  featsizes(toolarge) = [];
		end
		if (max(size(learnsizes)) > 1)
			error('For a feature curve, specify a scalar LEARNSIZE.');
		end;
	  featsizes = featsizes(:)';
		sizes = featsizes;
		featcurve = 1;
	end;

	% Fill the error structure.

  nw = length(classf(:));
  datname = getname(a);

  e.error   = zeros(nw,length(sizes));
  e.std     = zeros(nw,length(sizes));
  e.xvalues = sizes(:)';
  e.n       = nreps;
  e.names   = [];
  if (nreps > 1)
  	e.ylabel= ['Averaged error (' num2str(nreps) ' experiments)'];
  elseif (nreps == 1)
  	e.ylabel = 'Error';
  else
		error('Number of repetitions NREPS should be >= 1.');
	end;

  if (~isempty(datname))
		if (isempty(t))
			 	e.title = ['Bootstrapped learning curve on ' datname];
		else
			 	e.title = ['Learning curve on ' datname];
	 	end
	end;

 	if (learncurve)
   	e.xlabel = 'Training set size';
 	else
 		e.xlabel = 'Feature size';
 	end;

  if (sizes(end)/sizes(1) > 20)
  	e.plot = 'semilogx'; 				% If range too large, use a log-plot for X.
  end
  
	% Report progress.
  
	s1 = sprintf('clevals: %i classifiers: ',nw);
	prwaitbar(nw,s1);

	% Store the seed, to reset the random generator later for different
	% classifiers.

  seed = rand('state');

	% Loop over all classifiers (with index WI).

  for wi = 1:nw

		% Assert that CLASSF{WI} is an untrained mapping.
  	isuntrained(classf{wi});

  	name = getname(classf{wi});
    e.names = char(e.names,name);
    prwaitbar(nw,wi,[s1 name]);

		% E1 will contain the error estimates.

  	e1 = zeros(nreps,length(sizes));

		% Take care that classifiers use same training set.

  	rand('state',seed); seed2 = seed;

		% For N repetitions...		
    s2 = sprintf('clevals: %i repetitions: ',nreps);
		prwaitbar(nreps,s2);

  	for i = 1:nreps

			prwaitbar(nreps,i,[s2 int2str(i)]);
			if (isempty(t))

  			% Bootstrap. Store the randomly permuted indices of samples of class
  			% CI to use in this training set in JR(CI,:).

    		JR = zeros(c,max(learnsizes));
    		for ci = 1:c
    			JC = findnlab(a,ci);
  				
  				% Necessary for reproducable training sets: set the seed and store
  				% it after generation, so that next time we will use the previous one.
    			rand('state',seed2);		

    			R = ceil(rand(1,max(learnsizes))*length(JC));
    			JR(ci,:) = JC(R)';

    			seed2 = rand('state'); 
    		end
	
				t = a;

			end

  		% Either the outer loop or the inner loop will be traversed just once,
			% depending on whether we want to find a learning curve or feature
			% size curve.

  		ii = 0;
			nlearns = length(learnsizes);
			s3 = sprintf('cleval: %i sizes: ',nlearns);
			prwaitbar(nlearns,s3);
      
    	for li = 1:nlearns
				nj = learnsizes(li);
        prwaitbar(nlearns,li,[s3 int2str(li) ' (' int2str(nj) ')']);
			  nfeatsizes = length(featsizes);
			  s4 = sprintf('clevals: %i feature sizes: ',nfeatsizes);
			  prwaitbar(nfeatsizes,s4);
        
  			for fi = 1:nfeatsizes

          prwaitbar(nfeatsizes,fi,[s4 int2str(fi) ' (' int2str(featsizes(fi)) ')']);
  				ii = ii + 1;

    			% Ja will contain the indices for this training set, Jt for
					% the tuning set.

    			Ja = []; 
      		for ci = 1:c
      			Ja = [Ja;JR(ci,1:nj)']; 
      		end;
					
					Jt = setdiff(1:m,Ja);
    				
    			% Train classifier CLASSF{WI} on this training set and 
  				% calculate error on tuning set.

					if (learncurve)
	        	e1(i,ii) = testc(a, ...
														 a(Ja,1:featsizes(fi))*classf{wi});
					else
	        	e1(i,ii) = testc(t(Jt,1:featsizes(fi)), ...
														 a(Ja,1:featsizes(fi))*classf{wi});
					end;

    		end
        prwaitbar(0);
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
