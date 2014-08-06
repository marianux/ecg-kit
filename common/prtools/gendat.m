%GENDAT Random sampling of datasets for training and testing
% 
%  [A,B,IA,IB] = GENDAT(X,N,SEED)
%  [A,B,IA,IB] = X*GENDAT([],N,SEED)
%  [A,B,IA,IB] = X*GENDAT(N,SEED)
%  [A,B,IA,IB] = GENDAT(X,ALF,SEED)
%  [A,B,IA,IB] = X*GENDAT([],ALF,SEED)
%  [A,B,IA,IB] = X*GENDAT(ALF,SEED)
% 
% INPUT
%   X      Dataset.
%   N,ALF  Number/fraction of objects to be selected (def: bootstrapping).
%          Alternatively a vector of numbers of objects for each class.
%   SEED   A state of the random number generation according to RANDRESET
%
% OUTPUT
%   A,B    Datasets
%   IA,IB  Original indices from the dataset X
%
% DESCRIPTION
% Generation of N objects from dataset X. They are stored in dataset A,
% the remaining objects in dataset B. IA and IB are the indices of the
% objects selected from X for A and B. The random object generation follows
% the class prior probabilities. So is the prior probability of a class is
% PA, then in expectation PA*N objects are selected from that class. If N
% is large or if one of the classes has too few objects in A, the number of
% generated objects might be less than N.
% 
% If N is a vector of sizes, exactly N(i) objects are generated for class i.
% Classes are ordered as given by GETLABLIST(A).  
%
% If the function is called without specifying N, the data set X is
% bootstrapped and stored in A. Not selected samples are stored in B.
%
% ALF should be a scalar < 1. For each class a fraction ALF of the objects
% is selected for A and the not selected objects are stored in B.
%
% If X is a cell array of datasets the command is executed for each
% dataset separately. Results are stored in cell arrays. For each dataset
% the random seed is reset, resulting in aligned sets for the generated
% datasets if the sets in X were aligned.
% 
% EXAMPLES 
% See PREX_PLOTC.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, GENSUBSETS, RANDRESET

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function [A,B,IA,IB] = gendat(varargin)

  argin = shiftargin(varargin,'vector');
	argin = setdefaults(argin,[],[],[]);
  if mapping_task(argin,'definition')
    A = define_mapping(argin,'generator','Data sampling');
    return
  end
  
  % execution
  [X,N,seed] = deal(argin{:});
  if isempty(seed)
    seed = randreset;
  end
  randreset(seed);
	% If the input is a cell array of datasets, apply this procedure
  % to the individual datasets.
	if (iscell(X))
    [A,B,IA,IB] = X*feval(mfilename,[],N,seed);
    return
  end

	% When required, get the right number of objects from the given
	% fraction ALF.
  if isdouble(X), Xdouble = true; else, Xdouble = false; end
	if ~isdatafile(X), X = prdataset(X); end
	X = setlablist(X); % remove empty classes first
	[m,k,c] = getsize(X);
	% we need at least one class below:
	unlabeled = 0;
	if c==0, 
	   X=cdats(X,1); 
	   c=1;
		 unlabeled = 1; % we need to correct for labeling at the end
	end

	R = classsizes(X);
	if ~isempty(N) && length(N) ~= 1 && length(N) ~= c
		error('Data size should be scalar or a vector matching the number of classes')
	end
	if ~islabtype(X,'crisp') 
		if numel(N) > 1
			prwarning(1,'Specification of numbers of objects per class not possible for given label type')
			N = sum(N);
		end
		if N < 1, N = ceil(N*m); end
	end
	if ~isempty(N) && all(N < 1) && islabtype(X,'crisp')
		%DXD it should also be possible to have a fraction for each of the
		%classes, I think...
		if length(N)==1
			N = ceil(N*R);
		else
			N = ceil(N(:).*R(:));
		end
	end

	% Depending if N (or ALF) is given, the objects are created using
	% subsampling or bootstrapping.
	IA = [];
	if (nargin < 2) || (isempty(N))			% Bootstrap
		for i=1:c
			J = findnlab(X,i);
			K = ceil(rand(R(i),1)*R(i));
			IA = [IA; J(K)];
		end
	else																% Subsampling
		if ~islabtype(X,'crisp')
			K = randperm(m);
			if (N > m)
				%DXD: I would like to have just a warning:
				%error('More objects requested than available.')
				prwarning(4,'More objects requested than available.')
				%N = m;
        K = repmat(K,1,ceil(N/m));
			end
			IA = K(1:N);
		else
			%p = X.prior;   % avoid warning
			if isempty(X,'prior')
				p = classsizes(X);
				p = p/sum(p);
			else
				p = getprior(X);
			end
			%p = getprior(X);
			N = genclass(N,p);
			for i=1:c
				J = findnlab(X,i);
				K = randperm(R(i));
        if (N(i) > R(i))
					%DXD: I would like to have just a warning:
					%error('More objects requested than available.')
					prwarning(1,'More objects requested than available in class %d.',i)
					%N(i) = R(i);
        end
        K = repmat(K,1,ceil(N(i)/R(i)));
				IA = [IA; J(K(1:N(i)))];
			end
		end
	end

	% Finally, extract the datasets:
	IB = [1:m]';
	IB(IA) = [];
	if unlabeled 
		X = setlabels(X,[]); % reset unlabeling
	end
	A = X(IA,:);
	B = X(IB,:);
  if Xdouble
    A = +A; B = +B;
  end

return;
