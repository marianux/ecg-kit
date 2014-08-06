%VOTEC Voting combining classifier
% 
%  W = VOTEC(V)
%  W = V*VOTEC
% 
% INPUT
%  V   Set of classifiers
%
% OUTPUT
%  W   Voting combiner
%
% DESCRIPTION
% If V = [V1,V2,V3,...] is a stacked set of classifiers trained for the
% same classes, W is the voting combiner: it selects the class with the
% highest vote of the base classifiers. This might also be used as
% A*[V1,V2,V3]*VOTEC in which A is a dataset to be classified.
%
% The direct classifier outputs D = B*W for a test set B are posterior
% probability estimates. D(i,j) = (v+1) / (n+c), in which v is the number
% of votes object i receives for the j-th class. n is the total number of
% classifiers, c the total number of classes.
%
% The base classifiers may be combined in a stacked way (operating in the
% same feature space) by V = [V1,V2,V3, ... ] or in a parallel way
% (operating in different feature spaces) by V = [V1;V2;V3; ... ]
%
% EXAMPLES
% PREX_COMBINING
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, PRODC, MAXC, MINC,
% MEDIANC, MEANC, AVERAGEC, STACKED, PARALLEL

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands 

% $Id: votec.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function w = votec(p1)

		type = 'vote';               % define the operation processed by FIXEDCC.
	name = 'Voting combiner';    % define the name of the combiner.

	% this is the general procedure for all possible
	% calls of fixed combiners handled by FIXEDCC
	if nargin == 0
		w = prmapping('fixedcc','combiner',{[],type,name});
	else
		w = fixedcc(p1,[],type,name);
	end

	if isa(w,'prmapping')
		w = setname(w,name);
	end

	return
