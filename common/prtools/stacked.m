%STACKED Combining classifiers in the same feature space
%
%	  WC = STACKED(W1,W2,W3,  ....) or WC = [W1,W2,W3, ...]
%	  WC = STACKED({W1,W2,W3, ...}) or WC = [{W1,W2,W3, ...}]
%	  WC = STACKED(WC,W1,W2,  ....) or WC = [WC,W2,W3, ...]
%
% INPUT
%	  W1,W2,W3  Set of classifiers
%
% OUTPUT
%	  WC        Combined classifier
%
% DESCRIPTION
% The base classifiers (or mappings) W1, W2, W3, ... defined in the same
% feature space are combined in WC. This is a classifier defined for the
% same number of features as each of the base classifiers and with the
% combined set of outputs. So, for three two class classifiers defined for
% the classes 'c1' and 'c2', a dataset A is mapped by D = A*WC on the outputs
% 'c1','c2','c1','c2','c1','c2', which are the feature labels of D. Note that
% classification by LABELD(D) finds for each vector in D the feature label
% of the column with the maximum value. This is equivalent to using the
% maximum combiner MAXC,
%
% Other fixed combining rules like PRODC, MEANC, and VOTEC can be applied by
% D = A*WC*PRODC. A trained combiner like FISHERC has to be supplied with
% the appropriate training set by AC = A*WC; VC = AC*FISHERC. So the
% expression VC = A*WC*FISHERC yields a classifier, not a dataset as with
% fixed combining rules. This classifier operates in the intermediate
% feature space, the output space of the set of base classifiers. A new
% dataset B has to be mapped to this intermediate space first by BC = B*WC
% before it can be classified by D = BC*VC. As this is equivalent to D =
% B*WC*VC, the total trained combiner is WTC = WC*VC = WC*A*WC*FISHERC. To
% simplify this procedure PRTools executes the training of a combined
% classifier by WTC = A*(WC*FISHERC) as WTC = WC*A*WC*FISHERC.
%
% It is also possible to combine a set of untrained classifiers, e.g. WC =
% [LDC NMC KNNC([],1)]*CLASSC, in which CLASSC takes care that all outputs
% will be transformed to appropriate posterior probabilities. Training of
% all base classifiers is done by WC = A*WC. Again, this may be combined
% with training of a combiner by WTC = A*(WC*FISHERC).
%
% EXAMPLES
% PREX_COMBINING

% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, MAXC, MINC, MEANC,
% MEDIANC, PRODC, FISHERC, PARALLEL

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function w = stacked(varargin)

		% No arguments given: just return map information.

	if (nargin == 0) 
		w = prmapping(mfilename,'combiner');
		return
	end

	% Single argument: should be a mapping or cell array of mappings.

  if (nargin == 1)
    v = varargin{1};
    % If V is a single mapping, process it directly.
    if (~iscell(v))
      if (~isa(v,'prmapping'))
 	 			error('Mapping expected.')
      end
      w = prmapping('stacked',getmapping_type(v),{v},getlabels(v));
      w = set(w,'size',getsize(v));
    else
			% If V is a cell array of mappings, call this function recursively.
      if (size(v,1) ~= 1)
        error('Row of cells containing mappings expected')
 		 	end
 	 		w = feval(mfilename,v{:});
    end
    return
	end

	% Multiple arguments, all of which are mappings: combine them.

	if (~(isa(varargin{1},'prdataset') | (isa(varargin{1},'double') & ~isempty(varargin{1}))))
 
		% Get the first mapping.
    v1 = varargin{1};
    if (isempty(v1))
      start = 3; 
			v1 = varargin{2}; 
    else
      start = 2;
    end

    ismapping(v1);														% Assert V1 is a mapping.
    k = prod(getsize_in(v1)); 
		labels = getlabels(v1); 
		type = getmapping_type(v1);

		% If V1 is already a stacked mapping without output conversion,
		% unpack it to re-stack.
		
    if (~strcmp(getmapping_file(v1),mfilename)) | getout_conv(v1) > 1
      v = {v1};
    else
      v = getdata(v1);
    end
  	
		% Now stack the second to the last mapping onto the first.
		
    for j = start:nargin
      v2 = varargin{j};
      if (~strcmp(type,getmapping_type(v2)))
        error('All mappings should be of the same type.')
      end
      if (getsize(v2,1) ~= k)
        error('Mappings should have equal numbers of inputs.')
      end
      v = [v {v2}];
      if ischar(labels)
        labels = char(labels,getlabels(v2));
      elseif iscell(labels)
        labels = {labels{:} getlabels(v2)};
      else
        labels = [labels;getlabels(v2)];
      end
    end
    w = prmapping('stacked',type,v,labels,k);

  elseif isdatafile(varargin{1})
    
    % datafile * stacked mapping: store it as postprocessing
    
    w = addpostproc(varargin{1},stacked(varargin{2:end}));
    
  else

  	% The first argument is a dataset: apply the stacked mapping.

    a = varargin{1};
    v = varargin{2};
    if (~isa(v,'prmapping'))
      error('Mapping expected as second argument.')
    end

		if nargin==2 & isstacked(v)
			n = length(v.data);
		else
			v = varargin(2:end);
			n = nargin-1;
		end
		
		% Calculate W, the output of the stacked mapping on A.
    w = [];
    if n > 2
      s = sprintf('Concatenating %i mappings: ',n);
      prwaitbar(n,s);
    end
    for j = 1:n
      if n > 2, prwaitbar(n,j,[s int2str(j)]); end
      b = a*v{j};
      % If, for a mapping to 1D (e.g. a 2-class discriminant)
      % more than 1 output is returned, truncate.
      if (size(v{j},2) == 1)
        b = b(:,1);
      end
      w = [w b];					% Concatenate the outputs.
    end
    if n > 2
      prwaitbar(0); 
    end
    
  end

return
