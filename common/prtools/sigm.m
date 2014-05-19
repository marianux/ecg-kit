%SIGM Sigmoid map
% 
%   W = W*SIGM
%   B = A*SIGM
%   W = W*SIGM([],SCALE)
%   B = SIGM(A,SCALE)
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
% SEE ALSO
% DATASETS, MAPPINGS, CLASSC

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: sigm.m,v 1.3 2007/04/13 09:31:44 duin Exp $

function out = sigm (a,scale)

		if (nargin < 2)
		prwarning(3,'no scale supplied, assuming 1');
		scale = []; 
	end
	
	% Depending on the type of call, return a mapping or sigmoid-mapped data.

	if (nargin == 0) | (isempty(a))
		w = prmapping(mfilename,'fixed',scale);
		w = setname(w,'Sigmoidal Mapping');
		out = w;
	elseif (isempty(scale))
		out = 1./(1+exp(-a));
	else
		out = 1./(1+exp(-a/scale));
	end

	return
