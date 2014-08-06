% DCSC Dynamic Classifier Selection Combiner
% 
%   V = DCSC(A,W,K,TYPE)
%   V = A*DCSC([],W,K,TYPE)
%   V = A*(W*DCSC(K,TYPE)) 
%   D = B*V
%     
% INPUT
%   A      Dataset used for training base classifiers as well as combiner
%   B      Dataset used for testing (executing) the combiner
%   W      Set of trained or untrained base classifiers
%   K      Integer, number of neighbors
%   TYPE   'soft' (default) or 'crisp'
% 
% OUTPUT
%   V   Trained Dynamic Classifier Selection
%   D   Dataset with prob. products (over base classifiers) per class
%
% DESCRIPTION
% This dynamic classifier selection takes for every object to be
% classified the K nearest neighbors of an evaluation set (training set) A
% and determines which classifier performs best over this set of objects.
% If the base classifiers (STACKED or PARALLEL) are untrained, A is used to
% train them as well.
%
% The selection of the best classifier can be made in a soft or in a crisp
% way. If TYPE is 'soft' (default) classifier confidences are averaged (see
% CLASSC), otherwise the best classifier is selected by voting.
% 
% REFERENCE
% G. Giacinto and F. Roli, Methods for Dynamic Classifier Selection
% 10th Int. Conf. on Image Anal. and Proc., Venice, Italy (1999), 659-664.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, STACKED, NAIVEBC, CLASSC, TESTD, LABELD

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function OUT = dcsc(varargin)

name = 'Dynamic Classifier Selection';
DefaultNumNeigh = 20;
DefaultType = 'soft';

argin = shiftargin(varargin,{'prmapping','scalar'});
argin = setdefaults(argin,[],[],[],[]);
[par1,par2,par3,par4] = deal(argin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Empty Call  DCSC, DCSC([],K,TYPE)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isempty(par1)
	if isempty(par3), par3 = DefaultType; end
	if isempty(par2), par2 = DefaultNumNeigh; end
	% If there are no inputs, return an untrained mapping.
	% (PRTools transfers the latter into the first)
	OUT = prmapping(mfilename,'combiner',{par2,par3});
	OUT = setname(OUT,name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Storing Base Classifiers:  W*DCSC, DCSC(W,K,TYPE)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ismapping(par1)
	if isempty(par3), par3 = DefaultType; end
	if isempty(par2), par2 = DefaultNumNeigh; end
	% call like OUT = DCSC(W,k) or OUT = W*DCSC([],k)
	% store trained or untrained base classifiers, 
	% ready for later training of the combiner
	BaseClassf = par1;
	if ~isparallel(BaseClassf) & ~isstacked(BaseClassf)
		error('Parallel or stacked set of base classifiers expected');
	end
	if ~(isa(par2,'double') & length(par2)==1 & isint(par2))
		error('Number of neighbors should be integer scalar')
	end
	OUT = prmapping(mfilename,'untrained',{BaseClassf,par2,par3});
	OUT = setname(OUT,name);
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Training Combiner (and base classifiers if needed):  
%   A*(W*DCSC), A*(W*DCSC([],K,TYPE)), DCSC(A,W,K,TYPE)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif isdataset(par1) & ismapping(par2) & ...
		(isuntrained(par2) | (istrained(par2) & (isstacked(par2) | isparallel(par2))))
	% call like OUT = DCSC(TrainSet,W,k,type) or TrainSet*(W*DCSC([],k,type))
	% (PRTools transfers the latter into the first)
	% W is a set of trained or untrained set classifiers.
	if isempty(par4), par4 = DefaultType; end
	if isempty(par3), par3 = DefaultNumNeigh; end
	TrainSet = par1;
	BaseClassf = par2;
	islabtype(TrainSet,'crisp'); % allow crisp labels only
	isvaldfile(TrainSet,1,2); % at least one object per class, 2 classes
	TrainSet = setprior(TrainSet,getprior(TrainSet,0));  % avoid many warnings
	if ~(isstacked(BaseClassf) | isparallel(BaseClassf))
		if iscombiner(BaseClassf)
			error('No base classifiers found')
		end
		Data = getdata(BaseClassf);
		BaseClassf = Data.BaseClassf; % base classifiers were already stored
  end
  if isuntrained(BaseClassf)              % base classifiers untrained, so train them!
	  BaseClassf = TrainSet*BaseClassf;
		n = length(getdata(BaseClassf));
	else  % base classifiers are trained, just check label lists
		BaseClassifiers = getdata(BaseClassf);
		n = length(BaseClassifiers);
		for j=1:n
    	if ~isequal(getlabels(BaseClassifiers{j}),getlablist(TrainSet))
     		error('Training set and base classifiers should deal with same labels')
    	end
		end
  end
  Data.BaseClassf = BaseClassf;
	
  if ~isempty(par3) || ~isempty(par4) % overrules previously defined k (NumNeigh)
    Data.NumNeigh = par3;
  end
  if ~isempty(par4) % overrules previously defined type
		if strcmp(par4,'soft')
			Data.SoftVote = 1;
		else
			Data.SoftVote = 0;
		end
  end
  % Let us determine for every object in the trainingset how good it is
  % classified by every base classifier 
  [m,p,c] = getsize(TrainSet);
  %n = size(BaseClassf,2)/c;  % nr. of base classifiers
  ClassTrain1 = reshape(+(TrainSet*BaseClassf),m,c,n);
	ClassTrain2 = zeros(m,n);
	for j=1:n
		ct = ClassTrain1(:,:,j)*normm;
		ct = ct(sub2ind([m,c],[1:m]',[getnlab(TrainSet)]));
		ClassTrain2(:,j) = ct;
	end
	if ~Data.SoftVote % find crisp outcomes
		[cmax,L] = max(ClassTrain2,[],2);
		ClassTrain2 = zeros(m,n);
		ClassTrain2(sub2ind([m,n],[1:m]',L)) = ones(length(L),1);
	end
	Data.ClassTrain = ClassTrain2;
	
	Data.TrainSet = TrainSet;
	OUT = prmapping(mfilename,'trained',Data,getlablist(TrainSet),size(TrainSet,2),getsize(TrainSet,3));
	OUT = setname(OUT,name);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Evaluation of the trained combiner V:  B*V, DCSC(B,V,K)
%   TYPE cannot be changed anymore
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif isdataset(par1) & ismapping(par2) & istrained(par2)
	if isempty(par3), par3 = DefaultNumNeigh; end
	TestSet = par1;
	BaseClassf = getdata(par2,'BaseClassf');
	TrainSet = getdata(par2,'TrainSet');
  [m,p,c] = getsize(TrainSet);
  n = size(BaseClassf,2)/c;  % nr. of base classifiers
  if ~isempty(par3) || ~isempty(par4)
    k = par3;
  else
	  k = getdata(par2,'NumNeigh');
  end
	
  D = distm(TestSet,TrainSet);  % disputable procedure for parallel base classifiers
  [dd,J] = sort(+D,2); % J stores the numeric labels of the nearest training objects
  ClassPerf = zeros(size(TestSet,1),n); % base classifier performances per testobject
  TestSize = size(TestSet,1);
  N = [1:TestSize];
  kk = k;
	ClassTrain = getdata(par2,'ClassTrain');
  while ~isempty(N) & kk > 0
    for j=1:n % find for every testobject and for all classifiers the mean 
      % performance (i.e. the mean of the correct class assignments) over 
      % its neighborhood
      ClassPerf(:,j) = mean(reshape(ClassTrain(J(:,1:kk),j),TestSize,kk),2); 
		end
    [cc,L(N,:)] = sort(-ClassPerf(N,:),2);
    NN = find(cc(:,1) == cc(:,2)); % solve ties
    N = N(NN);  % select objects that suffer from ties
    kk = kk-1;  % try once more with smaller set of neighbors
  end
  U = getdata(BaseClassf);
  d = zeros(TestSize,c);
	feats = 0;  % counter to find feature numbers for parallel classifiers
  for j=1:n
		featsize = size(U{j},1);
    Lj = find(L(:,1)==j);
		if ~isempty(Lj)
			if isparallel(BaseClassf) % retrieve features for parallel classifiers
				d(Lj,:) = +(TestSet(Lj,feats+1:feats+featsize)*U{j});
			else
				d(Lj,:) = +(TestSet(Lj,:)*U{j});
			end
		end
		feats = feats+featsize;
  end
  OUT = setdat(TestSet,d,par2);
  
else
  
  error('Illegal input');
  
end
  