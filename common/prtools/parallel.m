%PARALLEL Combining classifiers in different feature spaces
%
%    WC = PARALLEL(W1,W2,W3,  ....) or WC = [W1;W2;W3; ...]
%    WC = PARALLEL({W1;W2;W3; ...}) or WC = [{W1;W2;W3; ...}]
%    WC = PARALLEL(WC,W1,W2,  ....) or WC = [WC;W2;W3; ...]
%    WC = PARALELL(C);
%    WC = PARALLEL(WC,N);
%
% INPUT
%   W1,W2,...  Base classifiers to be combined.
%   WC         Parallel combined classifier
%   C          Cell array of classifiers
%   N          Integer array
%
% OUTPUT
%    WC        Combined classifier.
%
% DESCRIPTION
% The base classifiers (or mappings) W1, W2, W3, ... defined in different
% feature spaces are combined in WC. This is a classifier defined for the 
% total number of features and with the combined set of outputs. So, for three 
% two-class classifiers defined for the classes 'c1' and 'c2', a dataset A is 
% mapped by D = A*WC on the outputs 'c1','c2','c1','c2','c1','c2' which are 
% the feature labels of D. Note that classification by LABELD(D) finds for 
% each vector in D the feature label of the column with the maximum value.
% This is equivalent to using the maximum combiner MAXC.
%
% Other fixed combining rules like PRODC, MEANC, and VOTEC can be applied 
% by D = A*WC*PRODC etc. A trained combiner like FISHERC has to be supplied
% with the appropriate training set by AC = A*WC; VC = AC*FISHERC. So the
% expression VC = A*WC*FISHERC yields a classifier and not a dataset as with 
% fixed combining rules. This classifier operates in the intermediate feature 
% space, the output space of the set of base classifiers. A new dataset B has 
% to be mapped to this intermediate space first by BC = B*WC before it can be 
% classified by D = BC*VC. As this is equivalent to D = B*WC*VC, the total 
% trained combiner is WTC = WC*VC = WC*A*WC*FISHERC. To simplify this procedure
% PRTools executes the training of a combined classifier by 
% WTC = A*(WC*FISHERC) as WTC = WC*A*WC*FISHERC.
%
% In order to allow for training an untrained parallel combined classifier by
% A*WC the subsets of the features of A that apply for the individual base
% classifiers of WC should be known to WC. Subset sizes are stored in the 
% dataset description if constructed by horizontal concatenation: 
% A = [A1 A2 A3}; This can be overruled by storing them in the classifier
% by the call WC = PARALLEL(WC,N), in which N is an array of integers, such
% that sum(N) equals the feature size of A.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, MAXC, MINC, MEANC, MEDIANC, PRODC, FISHERC, STACKED
%
% EXAMPLES
% See PREX_COMBINING.

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: parallel.m,v 1.4 2009/03/18 16:17:59 duin Exp $

function out = parallel (varargin)

		% If there are no arguments, just return an empty map.
	if (nargin == 0) 
		out = prmapping(mfilename,'combiner');
		return
	end

	% If there is one argument, assume it is a mapping or cell array of these.
	if (nargin == 1)
		v = varargin{1};
		if (~iscell(v))
			ismapping(v);         % Assert that V is a single mapping.
			out = prmapping('parallel',getmapping_type(v),{v},getlabels(v));
			out = set(out,'size',getsize(v));
		else
			% If V is a cell array of mappings, unpack it and call this function
			% for each cell.
			out = feval(mfilename,v{:});
		end
		return
	end

	if ( nargin == 2 & ( isparallel(varargin{1}) | iscell(varargin{1}) ) ...
			& ~ismapping(varargin{2}) )
		% special case, store dataset sizes
		if iscell(varargin{1})
			w = parallel(varargin{1});
		else
			w = varargin{1};
		end
		n = varargin{2};
		if length(w.data) ~= length(n)
			error('Wrong number of classsizes')
		end
		w.data = [w.data {n}];
		out = w;
		
	% If there are multiple arguments and the first is not a dataset,
	% combine the supplied mappings.
	% elseif ((nargin > 2) | ~isa(varargin{1},'prdataset')) 
	elseif (~(isa(varargin{1},'prdataset')))

		v1 = varargin{1}; 
  	if (isempty(v1) & ~ismapping(v1))
  		start = 3; v1 = varargin{2}; 
  	else
  		start = 2;
  	end
    
		ismapping(v1);     % Assert that V1 is a mapping.
		[k,n]      = size(v1);             % Extract V1's characteristics.
		labels = getlabels(v1);
		type   = getmapping_type(v1);

		if (~strcmp(getmapping_file(v1),mfilename))
			% If V1 is not already a parallel combiner, make it into a cell array
			% of mappings.
			v = {v1};
		else		
			% V1 is already a parallel combiner: get the mapping data.
			v = getdata(v1);
		end
		
		% Now concatenate all base classifiers as cells in V.		
		for j = start:nargin
			v2 = varargin{j}; 
			if j==nargin & ~ismapping(v2)
				v = [v {v2}];
			else
				if (~strcmp(type,getmapping_type(v2)))
					error('mappings should be of equal type')
				end
				v      = [v {v2}]; 
				k      = k + size(v2,1);
				n      = n + size(v2,2);
				labels = [labels; getlabels(v2)];
			end
		end

		if length(v) == 1
			out = v{1};                                % just one mapping left: return it
		else
			out = prmapping('parallel',type,v,labels,k,n); % Construct the combined mapping.
		end
		
	else 

		% Execution: dataset * parallel_mapping.
		a = varargin{1}; 
		if nargin > 2 & ~ismapping(varargin{end})
			v = parallel(varargin(2:end-1),varargin{end});
		elseif nargin > 2
			v = parallel(varargin(2:end));
		else
			v = varargin{2}; 
		end
		ismapping(v);    % Assert that V is a mapping.
		out = []; n = 0;
		if isuntrained(v)
      % get feature sizes of datasets to be supplied to the classifiers
      s = v.data{end};
      if ismapping(s) 
        s = getuser(a,'featsets');
        if isempty(s)
          error(['Training of parallel combined untrained classifier not possible.' ...
						newline 'Feature sizes should be stored in the classifier first.'])
        end
      else
        v.data = v.data(1:end-1);
      end
 			if sum(s) ~= size(a,2)
				error('Feature size of dataset does not match with classifier')
      end
			r = length(s);
      if length(v.data) ~= r
        error('Number of classifiers does not match number of feature subsets')
      end
			t = sprintf('Training %i base classifiers: ',r);
			prwaitbar(r,t)
			for j=1:r
				prwaitbar(r,j,[t getname(v{j})]);
				N  = [n+1:n+s(j)];       % to features indexed by N.
				n = n + s(j);
				w  = a(:,N)*v{j};
				out = [out; w];
			end
			prwaitbar(0)
		else
			for j = 1:length(v.data)
				sz = size(v{j},1);     % Classifier V{J} is applied
				N  = [n+1:n+sz];       % to features indexed by N.
				b  = a(:,N)*v{j}; 
				if (size(v{j},2) == 1) % Restore 2D outputs for k->1 classifiers.
					b = b(:,1);
				end
				n = n + sz;
				out = [out b];         % Concatenate mapped data.
			end
		end
	end

	return
