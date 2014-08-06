%FEATSELLR Plus-L-takeaway-R feature selection for classification
% 
%  [W,RES] = FEATSELLR(A,CRIT,K,L,R,T)
%
% INPUT	
%   A     Dataset
%   CRIT  String name of the criterion or untrained mapping 
%         (optional; default: 'NN', i.e. 1-Nearest Neighbor error)
%   K     Number of features to select
%         (optional; default: return optimally ordered set of all features)
%   L     Number of features to select at a time (plus-L, default: 1), L ~= R
%   R     Number of features to deselect at a time (takeaway-R, default: 0)
%   T     Tuning set (optional)
%   N     Number of cross-validations (optional)
%
% OUTPUT
%   W     Output feature selection mapping
%   RES   Matrix with step-by-step results of the selection
%
% DESCRIPTION
% Floating selection of K features using the dataset A, by iteratively 
% selecting L optimal features and deselecting R. Starts from the full
% set of features when L < R, otherwise from the empty set. CRIT sets the 
% criterion used by the feature evaluation routine FEATEVAL. If the dataset
% T is given, it is used as a tuning set for FEATEVAL. Alternatively
% a number of cross-validations N may be supplied. For K = 0, the optimal
% feature set (maximum value of FEATEVAL) is returned. The result W can be
% used for selecting features by B*W. In this case, features are ranked
% optimally. 
% The selected features are stored in W.DATA and can be found by +W.
% In R, the search is reported step by step as:
% 
% 	RES(:,1)            : number of features
% 	RES(:,2)            : criterion value
% 	RES(:,3:3+max(L,R)) : added / deleted features
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, FEATEVAL, FEATSEL
% FEATSELO, FEATSELB, FEATSELF, FEATSELI, FEATSELP, FEATSELM, PRPROGRESS

% Copyright: D. de Ridder, D.deRidder@ewi.tudelft.nl
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

% $Id: featsellr.m,v 1.6 2009/11/27 08:53:00 duin Exp $

function [w,res] = featsellr(a,crit,ksel,l,r,t,fid)

		
	if (nargin < 2) | isempty(crit)
		prwarning(2,'No criterion specified, assuming 1-NN.');
		crit = 'NN';
	end
	if (nargin < 3) | isempty(ksel)
		ksel = 0; 		% Consider all the features and sort them.
	end
	if (nargin < 4) | isempty(l) 
		l = 1;
	end;
	if (nargin < 5) | isempty(r)
		r = 0;
	end;
	if (nargin < 6)
		prwarning(3,'No tuning set supplied.');
		t = [];
	end
	if (nargin < 7)
		fid = [];
	end

	if (l==r)
		error('In ''Plus-L-takeaway-R feature selection'' L should be unequal to R')
	end
	
	% No inputs arguments provided, return an untrained mapping.
	if (nargin == 0) | (isempty(a))
		w = prmapping('featsellr',{crit,ksel,t});
		w = setname(w,'+L-R FeatSel');
		return
	end

	isvaldfile(a,1,2);               % At least 1 object per class, 2 classes.
	a = testdatasize(a);
	a = setprior(a,getprior(a));
	if ~is_scalar(t), iscomdset(a,t); end

	[m,k,c] = getsize(a);
	featlist = getfeatlab(a);

    % Initialise counters.
	state.l = l; state.r = r;
	if (state.l > state.r)                  % Start from empty set.
		state.I = [];   			        % Pool of selected feature indices.
		state.critval_opt = 0;              % Maximum criterion value found so far.
		state.res = zeros(0,2+max(l,r));    % Report of selection process.
	else
   	state.I = [1:k];			                            % Pool of the feature indices to be used.
		state.critval_opt = feateval(a,crit,t);                 % Max. criterion found so far.
   	state.res = [k,state.critval_opt,zeros(1,max(l,r))];    % Report of selection process. 
	end;
	state.I_opt = state.I; state.critval = [];

    % Evaluate the criterion function for the entire dataset A.

	prprogress(fid,['\nfeatsellr: Plus-' num2str(l) '-takeaway-' num2str(r) ...
			' feature selection\n'])
	if ksel == 0
		if (l > r) 
			s = sprintf('Forward search of optimal feature set out of %i: ',k);
		else
			s = sprintf('Backward search of optimal feature set out of %i: ',k);
		end
		nsize = k;
	else
		s = sprintf('Selection of %i features: ',ksel);
		nsize = abs(length(state.I)-ksel);
	end
	prwaitbar(nsize,s);
	
	while (ksel == 0 & (((state.l > state.r) & (length(state.I) <= k-state.l)) | ...
                        ((state.l < state.r) & (length(state.I) > state.r)))) | ... 
          (ksel > 0  & state.l > state.r & length(state.I) < ksel) | ...
          (ksel > 0  & state.l < state.r & length(state.I) > ksel)

        % Calculate criterion values C for the features in I and find the feature 
        % giving the maximum. 
		if (l > r)
			state = plusl(a,crit,t,state,fid);
			if (r > 0), state = takeawayr(a,crit,t,state,fid); end;
			prwaitbar(nsize,length(state.I),[s int2str(length(state.I))]);
		else
			state = takeawayr(a,crit,t,state,fid);
			if (l > 0), state = plusl(a,crit,t,state,fid); end;
			prwaitbar(nsize,k-length(state.I),[s int2str(length(state.I))]);
		end;    
        
	end;
	prwaitbar(0);
    
    % Make sure we end up with exactly KSEL features.

	if (ksel == 0) & (state.l > state.r) & (length(state.I) < k)
		state.l = k - length(state.I);
		state = plusl(a,crit,t,state,fid);
	elseif (ksel == 0) & (state.l < state.r) & (length(state.I) > 1)
		state.r = length(state.I) - 1;
		state = takeawayr(a,crit,t,state,fid);
	elseif (ksel > 0) & (state.l > state.r) & (length(state.I) > ksel)
		state.r = length(state.I) - ksel;
		state = takeawayr(a,crit,t,state,fid);
	elseif (ksel > 0) & (state.l < state.r) & (length(state.I) < ksel)
		state.l = ksel - length(state.I);
		state = plusl(a,crit,t,state,fid);
	end;
    
	if (ksel ~= 0), state.I_opt = state.I; end;

	w = featsel(k,state.I_opt);
  w = setmapping_type(w,'trained');
  w = setsize(w,[k length(state.I_opt)]);
	if ~isempty(featlist)
		w = setlabels(w,featlist(state.I_opt,:));
	end
	w = setname(w,'+L-R FeatSel');

	res = state.res;
	prprogress(fid,'featsellr finished\n')
	
return;

function state = plusl(a,crit,t,state,fid)

    % J contains the indices of the nonselected features.
	[m,k,c] = getsize(a); J = setdiff(1:k,state.I);

	critval_opt = 0; Jsub_opt = [];

    % Find all possible choices of L indices out of J.    
	[Jsub,ind] = npickk(J,state.l);
	while (~isempty(Jsub))
        % Calculate the criterion when subset JSUB is added.
		if (isempty(t))
			critval = feateval(a(:,[state.I Jsub]),crit);
		elseif is_scalar(t)
			critval = feateval(a(:,[state.I Jsub]),crit,t);
		else
			critval = feateval(a(:,[state.I Jsub]),crit,t(:,[state.I Jsub]));
		end;
        % Store the best subset and its criterion value.
		if (critval > critval_opt)
			critval_opt = critval;
			Jsub_opt    = Jsub;
		end;
		[Jsub,ind] = npickk(J,state.l,ind);
	end;
    
	state.I = [state.I Jsub_opt];

	if (critval_opt > state.critval_opt) | ...
			(critval_opt == state.critval_opt & ...
			length(state.I) < length(state.I_opt)), 
		state.critval_opt = critval_opt; 
		state.I_opt = state.I;  
	end;
    
	line = [length(state.I),critval_opt,...
				[Jsub_opt zeros(1,size(state.res,2)-2-length(Jsub_opt))]];
	prprogress(fid,'  %d %f \n',line(1:2)); 
	state.res = [state.res; line];

return

function state = takeawayr(a,crit,t,state,fid)

    % J contains the indices of the selected features.
	[m,k,c] = getsize(a); J = state.I;

	critval_opt = 0; Jsub_opt = [];

    % Find all possible choices of L indices out of J.    
	[Jsub,ind] = npickk(J,state.r);
	while (~isempty(Jsub))
        % Calculate the criterion when subset JSUB is removed.
		if (isempty(t))
			critval = feateval(a(:,[setdiff(state.I,Jsub)]),crit);
		elseif is_scalar(t)
			critval = feateval(a(:,[setdiff(state.I,Jsub)]),crit,t);
		else
			critval = feateval(a(:,[setdiff(state.I,Jsub)]),crit,...
                               t(:,[setdiff(state.I,Jsub)]));
		end;
        % Store the best subset and its criterion value.
		if (critval > critval_opt)
			critval_opt = critval;
			Jsub_opt    = Jsub;
		end;
		[Jsub,ind] = npickk(J,state.r,ind);
	end;
    
	state.I = setdiff(state.I,Jsub_opt);

	if (critval_opt > state.critval_opt) | ...
			(critval_opt == state.critval_opt & ...
			length(state.I) < length(state.I_opt)), 
		state.critval_opt = critval_opt; 
		state.I_opt = state.I;  
	end;

	line = [length(state.I),critval_opt,...
            [-Jsub_opt zeros(1,size(state.res,2)-2-length(Jsub_opt))]];
	prprogress(fid,'  %d %f \n',line(1:2)); 
	state.res = [state.res; line];
    
return

% NPICKK - A kind-of-recursive implementation of NCHOOSEK. 
%
% Picks K elements out of vector V and returns them in R, with their indices
% in I. Initialise with [R,I] = NPICKK(V,K); subsequent calls should be
% [R,I] = NPICKK(V,K,I). If the set is exhausted, returns R = I = [].

function [r,i] = npickk(v,k,i)

    n = length(v);

    if (nargin < 3) | (isempty(i))
        i = 1:k; j = 1;
    else
        % Increase the last index.
        j = k; i(j) = i(j) + 1; 
        % While there's an overflow, increase the previous index
        % and reset the subsequent ones.
        while (i(j) > n-(k-j)) & (j >= 2)
            i(j-1) = i(j-1) + 1; 
            for l = j:k
                i(l)   = i(l-1) + 1;
            end;
            j = j - 1;
        end;
    end;

    % If the last index is too large, we're done.
    if (i(end) > n)
        r = []; i = [];
    else
        r = v(i);
    end;
    
return

