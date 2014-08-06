%FEATSELI Trainable mapping for individual feature selection
% 
%   [W,R] = FEATSELI(A,CRIT,K,T)
%   [W,R] = A*FEATSELI([],CRIT,K,T)
%   [W,R] = A*FEATSELI(CRIT,K,T)
%   [W,R] = FEATSELI(A,CRIT,K,N)
%   [W,R] = A*FEATSELI([],CRIT,K,N)
%   [W,R] = A*FEATSELI(CRIT,K,N)
% 
% INPUT	
%   A    Training dataset
%   CRIT Name of the criterion or untrained mapping
%        (default: 'NN', i.e. the LOO 1-Nearest Neighbor error)
%   K    Number of features to select (default: sort all features)
%   T    Tuning dataset (optional)
%   N    Number of cross-validations (optional)
%
% OUTPUT
%   W    Feature selection mapping
%   R    Matrix with criterion values
%
% DESCRIPTION
% Individual selection of K features using the dataset A. CRIT sets the
% criterion used by the feature evaluation routine FEATEVAL. If the dataset
% T is given, it is used as test set for FEATEVAL. For K = 0 all features are
% selected, but reordered according to the criterion. The result W can be
% used for selecting features using B*W. 
% The selected features are stored in W.DATA and can be found by +W.
% In R, the search is reported step by step as:
% 
% 	R(:,1) : number of features
% 	R(:,2) : criterion value
% 	R(:,3) : added / deleted feature
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, FEATEVAL, FEATSELO, FEATSELB, FEATSELF,
% FEATSEL, FEATSELP,FEATSELM

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: featseli.m,v 1.5 2009/07/01 07:48:05 davidt Exp $

function [w,r] = featseli(varargin)

  varargin = shiftargin(varargin,{'char','prmapping'});
  argin = setdefaults(varargin,[],'NN',0,[]);
  if mapping_task(argin,'definition')
    w = define_mapping(argin,'untrained','Individual FeatSel');
    return
  end
    
  [a,crit,ksel,t] = deal(argin{:});

	[m,k,c] = getsize(a); featlist = getfeatlab(a);

	% If KSEL is not given, return all features.

	if (ksel == 0), ksel = k; end
	
	isvaldfile(a,1,2); % at least 1 object per class, 2 classes
	a = testdatasize(a);
	if isdataset(t), iscomdset(a,t); end
	
	critval = zeros(k,1);
	
	% Evaluate each feature in turn.

	s = sprintf('Evaluation of %i features: ',k);
	prwaitbar(k,s);
  a = setprior(a,getprior(a));
	if (isempty(t))
		for j = 1:k
			prwaitbar(k,j,[s int2str(j)]);
			critval(j) = feateval(a(:,j),crit);
		end
  elseif isdataset(t)
		for j = 1:k
			prwaitbar(k,j,[s int2str(j)]);
			critval(j) = feateval(a(:,j),crit,t(:,j));
    end
  else
		for j = 1:k
			prwaitbar(k,j,[s int2str(j)]);
			critval(j) = feateval(a(:,j),crit,t);
    end
	end
	prwaitbar(0);

	% Sort the features by criterion value (maximum first).
	[critval_sorted,J] = sort(-critval);
	r = [[1:k]', -critval_sorted, J];
	J = J(1:ksel)'; 

	% Return the mapping found.
	w = featsel(k,J);
  w = setmapping_type(w,'trained');
  w = setsize(w,[k length(J)]);
	if ~isempty(featlist)
		w = setlabels(w,featlist(J,:));
	end
	w = setname(w,'Individual FeatSel');

return
