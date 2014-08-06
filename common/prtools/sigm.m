%SIGM Sigmoid map
% 
%   W = W*SIGM(SCALE)
%
%   B = SIGM(A,SCALE)
%   B = A*SIGM([],SCALE)
%   B = A*SIGM(SCALE)
% 
% INPUT
%   A        Dataset (optional)
%   SCALE    Scaling parameter (optional, default: 1)
% 
% OUTPUT
%   W        Sigmoid mapping, or
%   B        Dataset A mapped by sigmoid mapping
%
% DESCRIPTION
% Sigmoidal transformation, useful to transform a map to classifier,
% producing posterior probability estimates. The parameter SCALE scales the
% data first (A/SCALE), before the transformation. Default: SCALE = 1, i.e.
% no scaling.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, CLASSC

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

function out = sigm (varargin)

  argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],1);
  
  if mapping_task(argin,'definition')
    out = define_mapping(argin,'fixed','Sigmoidal Mapping');
    
  else			% Evaluate
    
    [a,scale] = deal(argin{:});
		out = 1./(1+exp(-a/scale));
    
	end

	return
