%FEATSELO Branch and bound feature selection
% 
%   [W,R] = FEATSELO(A,CRIT,K,T)
%   [W,R] = A*FEATSELO([],CRIT,K,T)
%   [W,R] = A*FEATSELO(CRIT,K,T)
%   [W,R] = FEATSELO(A,CRIT,K,N)
%   [W,R] = A*FEATSELO([],CRIT,K,N)
%   [W,R] = A*FEATSELO(CRIT,K,N)
%
% INPUT	
%   A     Input dataset
%   CRIT  String name of the criterion or untrained mapping 
%           (optional, def= 'maha-s')
%   K     Numner of features to select (optional, def: K=2)
%   T     Validation set (optional)
%   N     Number of cross-validations (optional)
%
% OUTPUT
%   W     Output feature selection mapping
%   R     Matrix with step-by-step results
% 
% DESCRIPTION
% Backward selection of K features by baktracking using the branch 
% and bound procedure on the data set A. CRIT sets the criterion 
% used by the feature evaluation routine FEATEVAL. If the data set T 
% is given, it is used as test set for FEATEVAL. Alternatively a number
% of cross-validations N may be supplied. The resulting W can be used for
% the selecting features of a dataset B by B*W. 
% The selected features are stored in W.DATA and can be found by +W.
% 
% This procedure finds the optimum feature set if a monotoneous 
% criterion is used. The use of a testset does not guarantee that.
%
% REFERENCE
% P. M. Narendra and K. Fukunaga
% A Branch and Bound Algorithm for Feature Subset Selection,
% IEEE Trans. Computer, 26(9), pp. 917-922, September 1977
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, FEATEVAL, FEATSELF, FEATSELB, FEATSELI,
% FEATSEL, FEATSELP, FEATSELM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: featselo.m,v 1.7 2009/07/01 09:33:23 duin Exp $

function [W,R] = featselo(varargin)

  varargin = shiftargin(varargin,{'char','prmapping'});
  argin = setdefaults(varargin,[],'maha-s',2,[],[]);
  if mapping_task(argin,'definition')
    W = define_mapping(argin,'untrained','B&B FeatSel');
    return
  end
    
  [A,crit,kmin,T,fid] = deal(argin{:});

	isvaldfile(A,1,2); % at least 1 object per class, 2 classes
	A = testdatasize(A);
	if isdataset(T), iscomdset(A,T); end

	[m,k,c] = getsize(A);
	featlist = getfeatlab(A);
  A = setprior(A,getprior(A));

	if ((kmin < 1) | (kmin >= k))
		error('The desired feature size should be > 0 and < dataset feature size')
	end
	
	% space for criteria values
	feat = zeros(1,k);

   % Get performance of the individual features:
	if isempty(T)
		for j=1:k
			feat(j) = feateval(A(:,j),crit);
		end
	elseif is_scalar(T)
		for j=1:k
			feat(j) = feateval(A(:,j),crit,T);
		end
	else
		for j=1:k
			feat(j) = feateval(A(:,j),crit,T(:,j));
		end
	end

   % Get the kmin worst(?) individual features according to their
   % individual performance:
	[F,S] = sort(feat);
	
	%sometimes the above line is bad compared to the following two
	%w = featselb(A,crit,[]);
	%S = fliplr(+w);
	
	Iopt = [k-kmin+1:k];

	I = [1:k];
	J = [zeros(1,kmin),1:(k-kmin-1),k-kmin+1,k+1];
	level = k;

   % Get the performance of Iopt
	if isdataset(T)
		bound = feateval(A(:,S(Iopt)),crit,T(:,S(Iopt)));
	elseif is_scalar(T)
		bound = feateval(A(:,S(Iopt)),crit,T);
	else
		bound = feateval(A(:,S(Iopt)),crit);
	end

	C = inf;
	prwaitbar(100,'Branch & Bound Feature Selection')
	iter = 0;
	while numel(I) > 0 && J(k+1) == k+1;
		iter = iter+1; 
		prwaitbar(100,100-100*exp(-iter/25),['Branch & Bound Feature Selection: ' num2str(iter)]);
		if J(level) == J(level+1) | level <= kmin | C <= bound
			J(level) = level - kmin;
			level = level + 1;
			I = sort([I,J(level)]);
			J(level) = J(level) + 1;
			C = inf;
		else
			I(J(level)) = [];
			level = level - 1;
			if J(level+1) < 3 & level == kmin+1 & 0 % never happens ??
				;
			else
				if isdataset(T)
					C = feateval(A(:,S(I)),crit,T(:,S(I)));
				elseif is_scalar(T)
					C = feateval(A(:,S(I)),crit,T);
				else
					C = feateval(A(:,S(I)),crit);
        end
				if level == kmin & C > bound
					bound = C;
					Iopt = I;
          disp([bound,iter,numel(I)])
        end
			end
		end
  end
	prwaitbar(0);

   % Store the optimal features in the mapping:
	W = featsel(k,S(Iopt));
  W = setmapping_type(W,'trained');
  W = setsize(W,[k length(S(Iopt))]);
	if ~isempty(featlist)
		W = setlabels(W,featlist(S(Iopt),:));
	end
	W = setname(W,'B&B FeatSel');

	R = [];  %DXD I'm still not sure what to return

return
