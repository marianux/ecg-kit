%MCLASSM Computation of a combined, multi-class based mapping
%
%  W = MCLASSM(A,MAPPING,MODE,PAR)
%  W = A*MCLASSM([],MAPPING,MODE,PAR)
%  W = A*MCLASSM(MAPPING,MODE,PAR)
%
% INPUT
%   A       Dataset
%   MAPPING Untrained mapping	
%   MODE    Combining mode (optional; default: 'weight')
%   PAR     Parameter needed for the combining
%
% OUTPUT
%   W       Combined mapping
%
% DESCRIPTION
% If A is a unlabeled dataset or double matrix, it is converted to a 
% one-class dataset. For one-class datasets A, the mapping is computed,
% calling the untrained MAPPING using the labeled samples of A only. 
% For multi-class datasets separate mappings are determined for each class
% in A. They are combined as defined by MODE and PAR. The following
% combining rules are supported:
% 'weight': weight the mapping outcome for class j by PAR(j) and sum
%           over the classes. This is useful for densities in which case
%           PAR is typically the set of class priors (these are in fact
%           the defaults if MODE = 'weight').
% 'mean'  : combine by averaging.
% 'min'   : combine by the minimum rule.
% 'max'   : combine by the maximum rule.
%
% This routine is only defined for datasets with crisp labels.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function w = mclassm(varargin)

  argin = shiftargin(varargin,'prmapping');
  argin = setdefaults(argin,[],[],'weight',[]);
  
  if mapping_task(argin,'definition')
    w = define_mapping(argin,'untrained');
    return
  end
    
  [a,mapp,mode,par] = deal(argin{:});	
	
	if ~isa(mapp,'prmapping') | ~isuntrained(mapp)
		error('Second parameter should be untrained mapping')
	end

	[a,c,lablist,p] = cdats(a);
	
	%islabtype(a,'crisp');
	
	if c == 1
		w = a*mapp;
	else
		w = [];
		for j=1:c
			b = seldat(a,j);
			w = [w,b*mapp];
		end
		if ismapping(mode)
			w = w*mode;
		elseif isstr(mode)
			switch mode
			case 'mean'
				w = w*meanc;
			case 'min'
				w = w*minc;
			case 'max'
				w = w*maxc;
			case 'weight'
				if size(w{1},2) ~= 1
					error('Weighted combining only implemented for K to 1 mappings')
				end
				if isempty(par)
					par = p;
				end
				lenw = length(w.data);
				if length(par) ~= lenw
					error(['Wrong number of weights. It should be ' num2str(lenw)])
				end
				w = w*affine(par(:));
			otherwise
				error('Unknown combining mode')
			end
		end
		
		w = setsize(w,[size(a,2),1]);
		w = setlabels(w,lablist);
		
	end
	
return

