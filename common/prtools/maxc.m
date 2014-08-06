%MAXC Maximum combining classifier
% 
% 	W = MAXC(V)
% 	W = V*MAXC
% 
% INPUT
%   V   Stacked set of classifiers
%
% OUTPUT
%   W   Combined classifier using max-rule
%
% DESCRIPTION
% If V = [V1,V2,V3, ... ] is a set of classifiers trained on the same
% classes, then W is the maximum combiner: it selects the class that gives
% the maximal output of the input classifiers. This might also be used as
% A*[V1,V2,V3]*MAXC in which A is a dataset to be classified. Consequently,
% if S is a similarity matrix with class feature labels (e.g. S =
% A*PROXM(A,'r')), then S*MAXC*LABELD is the nearest neighbor classifier.
% 
% If it is desired to operate on posterior probabilities then the input
% classifiers should be extended to output these, using V = V*CLASSC.
%
% The base classifiers may be combined in a stacked way (operating in the
% same feature space by V = [V1,V2,V3, ... ] or in a parallel way (operating
% in different feature spaces) by V = [V1;V2;V3; ... ].
% 
% EXAMPLES
% see PREX_COMBINING
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, VOTEC, MINC, MEANC, MEDIANC, PRODC, AVERAGEC,
% STACKED, PARALLEL

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands 

% $Id: maxc.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function w = maxc(p1)

	% This is the general procedure for all possible calls of fixed combiners:
	%   they are handled by FIXEDCC.

	type = 'max';              	% Define the operation processed by FIXEDCC.
	name = 'Maximum combiner'; 	% Define the name of the combiner.

  % Possible calls: MAXC, MAXC(W) or MAXC(A,W).
                           
	if (nargin == 0)						% No arguments given: return untrained mapping.
		w = prmapping('fixedcc','combiner',{[],type,name});
	else
		w = fixedcc(p1,[],type,name);
  end

  if (isa(w,'prmapping'))
  	w = setname(w,name);
  end

return
