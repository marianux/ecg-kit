%FEATEVAL Evaluation of feature set for classification
% 
% 	J = FEATEVAL(A,CRIT,T)
%   J = A*FEATEVAL([],CRIT,T)
%   J = A*FEATEVAL(CRIT,T)
% 	J = FEATEVAL(A,CRIT,N)
%   J = A*FEATEVAL([],CRIT,N)
%   J = A*FEATEVAL(CRIT,N)
% 
% INPUT
%       A      input dataset
%       CRIT   string name of a method or untrained mapping, default 'NN'
%       T      validation dataset (optional)
%       N      number of cross-validations (optional)
%
% OUTPUT
%       J      scalar criterion value
%
% DESCRIPTION
% Evaluation of features by the criterion CRIT, using objects in the
% dataset A. The larger J, the better. Resulting J-values are
% incomparable over the various methods.
% The following methods are supported:
%  
%   crit='in-in' : inter-intra distance.
%   crit='maha-s': sum of estimated Mahalanobis distances.
%   crit='maha-m': minimum of estimated Mahalanobis distances.
%   crit='eucl-s': sum of squared Euclidean distances.
%   crit='eucl-m': minimum of squared Euclidean distances.
%   crit='NN'    : 1-Nearest Neighbour leave-one-out
%                  classification performance (default).
%                  (performance = 1 - error). 
%   crit='mad'   : mean absolute deviation (only for regression!)
%   crit='mse'   : mean squared error (only for regression!)
% 
% For classification problems, CRIT can also be any untrained
% classifier, e.g. LDC([],1e-6,1e-6). Then the classification error is
% used for a performance estimate. If supplied, the dataset T is used
% for obtaining an unbiased estimate of the performance of classifiers
%trained with the dataset A. If a number of cross-validations N is
% supplied, the routine is run for N times with different training and
% test sets generated from A by cross-validation. Results are averaged.
% If T nor N are given, the apparent performance on A is used. 
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, FEATSELO, FEATSELB, FEATSELF, FEATSELP, FEATSELM, FEATRANK

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% REVISIONS
% DXD1: David Tax, 08-05-2003
%       I added the inter/intra distance criterion.
% DXD2: David Tax, 10-06-2011
%       I added the MAD and MSE criteria for regression.
% RD  : Conversion to mapping, 08-08-2012

function J = feateval(varargin)
	
argin = shiftargin(varargin,'char');
argin = setdefaults(argin,[],'NN',[]);
if mapping_task(argin,'definition')
  J = define_mapping(argin,'fixed');
else
  [a,crit,t] = deal(argin{:});
	[ma,k,c] = getsize(a);

	if is_scalar(t) & ~isdataset(t) % cross-validation desired, t rotations
		if ~ismapping(crit) | ~isuntrained(crit)
 			error('Cross-validation only possible with untrained classifiers')
		end
		J = 1-prcrossval(a,crit,t);
		return
	end
		
%	islabtype(a,'crisp');
	isvaldfile(a,1,2); % at least 1 object per class, 2 classes
	a = testdatasize(a);
	iscomdset(a,t);

	if isstr(crit)
		%DXD1
		if strcmp(crit,'in-in')     % inter/intra distances
			islabtype(a,'crisp','soft');
			if isempty(t)
				[U,G] = meancov(a);
			else
				[U,G] = meancov(t);
			end
			S_b = prcov(+U); % between scatter
			prior = getprior(a);
			S_w = reshape(sum(reshape(G,k*k,c)*prior',2),k,k); % within scatter
			J = trace(prinv(S_w)*S_b);
		elseif strcmp(crit,'maha-s') | strcmp(crit,'maha-m') % Mahalanobis distances
			islabtype(a,'crisp','soft');
			if isempty(t)
				D = distmaha(a);
			else
				[U,G] = meancov(a);
				D = distmaha(t,U,G);
				D = meancov(D);
			end
			if strcmp(crit,'maha-m')
				D = D + realmax*eye(c);
				J = min(min(D));
			else
				J = sum(sum(D))/2; 
			end
		elseif strcmp(crit,'eucl-s') | strcmp(crit,'eucl-m') % Euclidean distances
			islabtype(a,'crisp','soft');
			U = meancov(a);
			if isempty(t)
				D = distm(U);
			else
				D = distm(t,U);
				D = meancov(D);
			end
			if strcmp(crit,'eucl-m')
				D = D + realmax*eye(c);
				J = min(min(D));
			else
				J = sum(sum(D))/2; 
			end
		elseif strcmp(crit,'NN')	% 1-NN performance
			islabtype(a,'crisp','soft');
			if isempty(t)
				J = 1 - testk(a,1);
			else
				J = 1 - testk(a,1,t);
			end
		elseif strcmp(crit,'kcentres') || strcmp(crit,'kcenters') % data radius, unsupervised
				% assumes disrep, so experimental
			J = max(min(+a,[],2));
			if J == 0, J = inf; else J = 1/J; end
		elseif strcmp(crit,'representation_error') || strcmp(crit,'reperror')% also unsupervised
			J = mean(min(+a,[],2));
			if J == 0, J = inf; else J = 1/J; end
		elseif strcmp(crit,'mad') % Mean Absolute Deviation for regression
			J = 1-testr(a*linearr(a,0.001,1),'mad');
		elseif strcmp(crit,'mse') % Mean Squared Error for regression
			J = 1-testr(a*linearr(a,0.001,1),'mse');
		else
			error('Criterion undefined');
		end
	else
		ismapping(crit);
		if isuntrained(crit)
			if isempty(t)
				J = 1 - (a * (a * crit) * testc);
			else
				J = 1 - (t * (a * crit) * testc);
			end
		elseif isfixed(crit)
			J = a * crit;
		else
			error('Criterion should be defined by an untrained or fixed mapping')
		end
	end
end
return
