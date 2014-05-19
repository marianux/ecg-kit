%FEATSELP Pudil's floating feature selection (forward)
% 
% [W,R] = FEATSELP(A,CRIT,K,T,FID)
%
% INPUT	
%   A    Training dataset
%   CRIT Name of the criterion or untrained mapping 
%        	(default: 'NN', 1-Nearest Neighbor error)
%   K    Number of features to select (default: K = 0, select optimal set)
%   T    Tuning dataset (optional)
%   N    Number of cross-validations (optional)
%   FID  File ID to write progress to (default [], see PRPROGRESS)
%
% OUTPUT
%   W    Feature selection mapping
%   R    Matrix with step-by-step results
% 
% DESCRIPTION
% Forward floating selection of K features using the dataset A. CRIT sets
% the criterion used by the feature evaluation routine FEATEVAL. If the
% dataset T is given, it is used as test set for FEATEVAL. Alternatively
% a number of cross-validations N may be supplied. For K = 0, the optimal
% feature set (maximum value of FEATEVAL) is returned. The result W can
% be used for selecting features in a dataset B using B*W.
% The selected features are stored in W.DATA and can be found by +W.
%
% Note: this routine is highly time consuming.
%
% In R the search is reported step by step:
% 
% 	R(:,1) : number of features
% 	R(:,2) : criterion value
% 	R(:,3) : added / deleted feature
% 
% SEE ALSO
% MAPPINGS, DATASETS, FEATEVAL, FEATSELO, FEATSELB, FEATSELI,
% FEATSEL, FEATSELF, FEATSELM, PRPROGRESS

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: featselp.m,v 1.5 2009/07/01 09:33:23 duin Exp $

function [w,r] = featselp(a,crit,ksel,t,fid)

		
	if (nargin < 2) | isempty(crit)
	prwarning(2,'no criterion specified, assuming NN');
		crit = 'NN'; 
	end
	if (nargin < 3) | isempty(ksel)
		ksel = 0; 
	end
	if (nargin < 4) | isempty(t)
	prwarning(3,'no tuning set supplied (risk of overfit)');
		t = [];
	end
	if (nargin < 5)
		fid = [];
	end

	% If no arguments are supplied, return an untrained mapping.

	if (nargin == 0) | (isempty(a))
		w = prmapping('featselp',{crit,ksel,t});
		w = setname(w,'Floating FeatSel');
		return
	end

	isvaldfile(a,1,2); % at least 1 object per class, 2 classes
	a = testdatasize(a);
	iscomdset(a,t);
	
	[m,k,c] = getsize(a); featlist = getfeatlab(a);

	% If KSEL is not given, return all features.

	if (ksel == 0)
		peak = 1; ksel = k; 
	else 
		peak = 0; 
	end

	if (~isempty(t))
		if (k ~= size(t,2))
			error('The sizes of the training and tuning dataset do not match.')
		end
	end

	critval_opt = zeros(1,k);	% Maximum criterion value for sets of all sizes.
	critval_max = 0;    			% Maximum criterion value found so far.
	I = [1:k];      			    % Pool of remaining feature indices.
	J = [];      			      	% Pool of selected feature indices.
	r = [];           			  % Result matrix with selection history.
	Iopt = J; 

	n = 0;
	prprogress(fid,'\nfeatselp: Pudils Floating Search\n')
	while (n < k)

		critval = zeros(1,length(I));

		% Add the best feature.

		for j = 1:length(I)
			L = [J,I(j)];					% Add one feature to the already selected ones.

			if (isempty(t))				% Evaluate the criterion function.
				critval(j) = feateval(a(:,L),crit);
			else
				critval(j) = feateval(a(:,L),crit,t(:,L));
			end

      % If this feature is the best so far and we have not yet selected
			% KSEL features, store it.

			if (critval(j) > critval_max) & (n < ksel)
				n_max = length(L);
				critval_max = critval(j);
				Iopt = L;
			end
		end

		[mx,j] = max(critval);  % Find best feature of the remaining ones,
		J = [J, I(j)];          %   add it to the set of selected features
		I(j) = [];              %   and remove it from the pool.

		% Store the best criterion value found for any set of n features.
		n = n + 1; critval_opt(n) = mx;

		r = [r; [n, mx, J(end)]];
		prprogress(fid,'  %d %f',r(end,1:2));
		prprogress(fid,' %i',J);
		prprogress(fid,'\n')
		
		% Now keep removing features until the criterion gets worse.

		while (n > 2)
			critval = zeros(1,n);
			for j = 1:n
				L = J; L(j) = [];		% Remove one feature from the selected ones.		

				if (isempty(t))			% Evaluate the criterion function.
					critval(j) = feateval(a(:,L),crit);
				else
					critval(j) = feateval(a(:,L),crit,t(:,L));
				end		

				% If removing this feature gives the best result so far (or
				% the same result using less features), and we have not yet
				% removed all KSEL features, store it.

				if ((critval(j) > critval_max) | ((critval(j) == critval_max) & ...
																					(length(L) < n_max))) & ...
					 (n <= ksel + 1)
					n_max = length(L);
					critval_max = critval(j);
					Iopt = L;
				end
			end

			% If this subset is better than any found before, store and report it.
			% Otherwise, stop removing features.

			[mx,j] = max(critval);

			if (mx > critval_opt(n-1))
				n = n - 1; critval_opt(n) = mx;
				I = [I,J(j)]; J(j) = [];
				r = [r; [n, mx, -I(end)]];
				prprogress(fid,'  %d %f',r(end,1:2));
				prprogress(fid,' %i',J);
				prprogress(fid,'\n')
			else
				break;
			end

		end

		% If we have found more than KSEL features, return the mapping using
		% the best KSEL features.

		if (n > ksel)
			if (ksel < length(Iopt))
				J = Iopt(1:ksel);
			else
				J = Iopt;
			end
			prprogress(fid,'featselp finished\n')
			w = featsel(k,J);
			if ~isempty(featlist)
				w = setlabels(w,featlist(J,:));
			end
			w = setname(w,'Floating FeatSel');
			return
		end

	end
	prprogress(fid,'featselp finished\n')
	
	% Return all features, sorted by their criterion value.

	w = featsel(k,Iopt);
	if ~isempty(featlist)
		w = setlabels(w,featlist(Iopt,:));
	end
	w = setname(w,'Floating FeatSel');

return
