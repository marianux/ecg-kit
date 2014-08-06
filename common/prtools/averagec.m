%AVERAGEC Combining of linear classifiers by averaging coefficients
% 
% 	W = AVERAGEC(V)
% 	W = V*AVERAGEC
%
% INPUT
% 	V  A set of affine base classifiers.
%
% OUTPUT
%  W  Combined classifier.
%
% DESCRIPTION
% Let V = [V1,V2,V3, ... ] is a set of affine classifiers trained on the same 
% classes, then W is the average combiner: it averages the coefficients of the 
% base classifiers, resulting in a new affine classifier. This might also be 
% used as A*[V1,V2,V3]*AVERAGEC, in which A is a dataset to be classified.
% 
% The base classifiers may be combined in a stacked way (operating in the same
% feature space by V = [V1,V2,V3, ... ] or in a parallel way (operating in 
% different feature spaces) by V = [V1;V2;V3; ... ].
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, VOTEC, MAXC, MINC, MEDIANC, PRODC, AVERAGEC, STACKED, 
% PARALLEL
% 
% EXAMPLES
% See PREX_COMBINING.

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands 

% $Id: averagec.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function w = averagec (p1)

  % The average combiner is constructed as a fixed combiner (FIXEDCC) of
  % type 'average'.

	type = 'average'; name = 'Average combiner';  

	% Possible calls: AVERAGEC, AVERAGEC(W) or AVERAGEC(A,W).

	if (nargin == 0)
		w = prmapping('fixedcc','combiner',{[],type,name});
	else
		w = fixedcc(p1,[],type,name);
	end

	if (isa(w,'prmapping'))
		w = setname(w,name);
	end

return
