%MINC Minimum combining classifier
% 
% 	W = MINC(V)
% 	W = V*MINC
% 
% INPUT
%   V     Set of classifiers
%
% OUTPUT
%   W    Minimum combining classifier on V
%
% DESCRIPTION
% If V = [V1,V2,V3, ... ] is a set of classifiers trained on the 
% same classes and W is the minimum combiner: it selects the class 
% with the minimum of the outputs of the input classifiers. This 
% might also be used as A*[V1,V2,V3]*MINC in which A is a dataset to 
% be classified. Consequently, if S is a dissimilarity matrix with
% class feature labels (e.g. S = A*PROXM(A,'d')) then S*MINC*LABELD
% is the nearest neighbor classifier.
% 
% If it is desired to operate on posterior probabilities then the 
% input classifiers should be extended like V = V*CLASSC;
%
% The base classifiers may be combined in a stacked way (operating
% in the same feature space by V = [V1,V2,V3, ... ] or in a parallel
% way (operating in different feature spaces) by V = [V1;V2;V3; ... ]
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, VOTEC, MAXC, MEANC, MEDIANC, PRODC,
% AVERAGEC, STACKED, PARALLEL
%
% EXAMPLES
% See PREX_COMBINING

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands 

% $Id: minc.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function w = minc(p1)

	type = 'min'; % define the operation processed by FIXEDCC.

	% define the name of the combiner. 
	% this is the general procedure for all possible calls of fixed combiners
	% handled by FIXEDCC
	name = 'Minimum combiner'; 

	if nargin == 0
		w = prmapping('fixedcc','combiner',{[],type,name});
	else
		w = fixedcc(p1,[],type,name);
	end

	if isa(w,'prmapping')
		w = setname(w,name);
	end

return
	
