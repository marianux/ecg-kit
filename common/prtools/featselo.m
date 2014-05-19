%FEATSELO Branch and bound feature selection
% 
% 	W = featselo(A,CRIT,K,T,FID)
%
% INPUT	
%   A     input dataset
%   CRIT  string name of the criterion or untrained mapping 
%           (optional, def= 'NN' 1-Nearest Neighbor error)
%   K     numner of features to select (optional, def: K=2)
%   T     validation set (optional)
%   N     Number of cross-validations (optional)
%   FID   File ID to write progress to (default [], see PRPROGRESS)
%
% OUTPUT
%   W     output feature selection mapping
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
% SEE ALSO
% MAPPINGS, DATASETS, FEATEVAL, FEATSELF, FEATSELB, FEATSELI,
% FEATSEL, FEATSELP, FEATSELM, PRPROGRESS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: featselo.m,v 1.7 2009/07/01 09:33:23 duin Exp $

function [W,R] = featselo(A,crit,kmin,T,fid)
		if nargin < 2 | isempty(crit), crit = 'NN'; end
	if nargin < 3 | isempty(kmin), kmin = 2; end
	if nargin < 4, T = []; end
	if (nargin < 5)
		fid = [];
	end

	if nargin == 0 | isempty(A)
      % Create an empty mapping:
		W = prmapping('featselo',{crit,kmin,T});
		W = setname(W,'B&B FeatSel');
		return
	end

	isvaldfile(A,1,2); % at least 1 object per class, 2 classes
	A = testdatasize(A);
	if ~is_scalar(T), iscomdset(A,T); end

	[m,k,c] = getsize(A);
	featlist = getfeatlab(A);

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
	prprogress(fid,'\nfeatselo : Branch & Bound Feature Selection\n')
	prwaitbar(100,'Branch & Bound Feature Selection')
	iter = 0;
	while length(I) > 0 & J(k+1) == k+1;
		iter = iter+1;
		prwaitbar(100,100-100*exp(-iter/1000));
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
				%prprogress(fid,'  level %i, crit: %10.5e',level,C);
				if level == kmin & C > bound
					bound = C;
					Iopt = I;
					prprogress(fid,'  level %i, crit: %10.5e',level,C);
					prprogress(fid,'  FeatSet: ');
					prprogress(fid,' %i',S(I));
					prprogress(fid,'\n')
				end
				%prprogress(fid,'\n')
			end
		end
	end
	prprogress(fid,'featselo  finished\n');
	prwaitbar(0);

   % Store the optimal features in the mapping:
	W = featsel(k,S(Iopt));
	if ~isempty(featlist)
		W = setlabels(W,featlist(S(Iopt),:));
	end
	W = setname(W,'B&B FeatSel');

	R = [];  %DXD I'm still not sure what to return

return
