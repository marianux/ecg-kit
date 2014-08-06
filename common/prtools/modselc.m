% MODSELC Model selection
% 
%   V = MODSELC(A,W,N,NREP)
%   V = A*(W*MODSELC([],N,NREP)) 
%     
% INPUT
%   A      Dataset used for training base classifiers and/or selection
%   B      Dataset used for testing (executing) the selector
%   W      Set of trained or untrained base classifiers
%   N      Number of crossvalidations, default 10
%   NREP   Number of crossvalidation repetitions
% 
% OUTPUT
%   V   Selected classidfer
%
% DESCRIPTION
% This routine selects out of a set of given classifiers stored in W the
% best one on the basis of N-fold crossvalidation (see CROSSVAL), which
% might be repeated NREP times. If W contains a set of already trained
% classifiers, N and NREP are neglected and just the best classifier
% according to the evaluation set A is returned.
%
% This routine can be considered as a classifier combiner based on global
% selection. See DCSC for local, dynamic selection.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, STACKED, DCSC, CLASSC, TESTD, LABELD

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function OUT = modselc(par1,par2,par3,par4)

name = 'Classifier Model Selection';
DefaultN = 10;
DefaultNrep = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Empty Call MODSELC, MODSELC([],N,NREP)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1 | isempty(par1)
	if nargin > 1 & (ismapping(par2) | iscell(par2))
		if nargin < 4 | isempty(par4), par4 = DefaultNrep; end
		if nargin < 3 | isempty(par3), par3 = DefaultN; end
		OUT = prmapping(mfilename,'untrained',{par2,par3,par4});
	else
		if nargin < 3 | isempty(par3), par3 = DefaultNrep; end
		if nargin < 2 | isempty(par2), par2 = DefaultN; end
		% If there are no inputs, return an untrained mapping.
		% (PRTools transfers the latter into the first)
		OUT = prmapping(mfilename,'combiner',{par2,par3});
	end
	OUT = setname(OUT,name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Storing Base Classifiers:  W*MODSELC, MODSELC(W,N,NREP)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ismapping(par1)
	if nargin < 3 | isempty(par3), par3 = DefaultNrep; end
	if nargin < 2 | isempty(par2), par2 = DefaultN; end
	% call like OUT = MODSELC(W,N) or OUT = W*MODSELC([],N)
	% store trained or untrained base classifiers, 
	% ready for later training of the combiner
	BaseClassf = par1;
	if ~isparallel(BaseClassf) & ~isstacked(BaseClassf)
		error('Parallel or stacked set of base classifiers expected');
	end
	OUT = prmapping(mfilename,'untrained',{BaseClassf,par2,par3});
	OUT = setname(OUT,name);
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Selection (and training base classifiers if needed):  
%   A*(W*MODSELC), A*(W*MODSELC([],N,NREP)), MODSELC(A,W,N,NREP)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif isdataset(par1) & ismapping(par2) & ...
		(isuntrained(par2) | (istrained(par2) & (isstacked(par2) | isparallel(par2))))
	% call like OUT = MODSELC(TrainSet,W,N,NREP) or TrainSet*(W*MODSELC([],N,NREP))
	% (PRTools transfers the latter into the first)
	% W is a set of trained or untrained set classifiers.
	if nargin < 4 | isempty(par4), par4 = DefaultNrep; end
	if nargin < 3 | isempty(par3), par3 = DefaultN; end
	TrainSet = par1;
	BaseClassf = par2;
	N = par3;
	Nrep = par4;
	isvaldfile(TrainSet,1,2); % at least one object per class, 2 classes
	TrainSet = setprior(TrainSet,getprior(TrainSet,0));  % avoid many warnings
	if ~(isstacked(BaseClassf) | isparallel(BaseClassf))
		if iscombiner(BaseClassf)
			error('No base classifiers found')
		end
		Data = getdata(BaseClassf);
		BaseClassf = Data.BaseClassf; % base classifiers were already stored
  end
  if isuntrained(BaseClassf)  % base classifiers untrained, crossval needed
		randstate = randreset; % makes routine reproducing
		if isparallel(BaseClassf)
			BaseClassf = getdata(BaseClassf);
			tsize = BaseClassf{end};
			if ismapping(tsize)
				error(['Training of parallel combined untrained classifier not possible.' ...
						newline 'Feature sizes should be stored in the classifier first.'])
			end
			csize = cumsum([0 tsize(:)']);
			n = length(BaseClassf)-1;
			e = zeros(1,n);
			rstate = randreset;
			s = sprintf('Crossvalidating %i base classifiers: ',n);
			prwaitbar(n,s)
			for j=1:n
				prwaitbar(n,j,[s getname(BaseClassf{j})]);
				randreset(1);
				e(j) = crossval(TrainSet(:,csize(j)+1:csize(j+1)),BaseClassf{j},N,Nrep);
			end
			prwaitbar(0);
			[ee,k] = min(e);
			J = [csize(k)+1:csize(k+1)];
			OUT = TrainSet(:,J)*BaseClassf{k};
			OUT = featsel(csize(end),J)*OUT;
		else
			if isstacked(BaseClassf)
				BaseClassf = getdata(BaseClassf);
			end
			randreset(1);
			e = crossval(TrainSet,BaseClassf,N,Nrep);
			[ee,k] = min(e);
			OUT = TrainSet*BaseClassf{k};
		end
		randreset(randstate);
	else  % base classifiers are trained, check label lists and find best one
		BaseClassf = getdata(BaseClassf);
		n = length(BaseClassf);
		for j=1:n
    	if ~isequal(getlabels(BaseClassf{j}),getlablist(TrainSet))
     		error('Training set and base classifiers should deal with same labels')
    	end
		end
		e = testc(TrainSet,BaseClassf);
		[ee,k] = min([e{:}]);
		OUT = BaseClassf{k};
  end
  
else
  
  error('Illegal input');
  
end
  