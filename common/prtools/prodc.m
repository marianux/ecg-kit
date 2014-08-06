%PRODC Product combining classifier
% 
%  W = PRODC(V)
%  W = V*PRODC
% 
% INPUT
%   V    Set of classifiers trained on the same classes
%
% OUTPUT
%   W    Product combiner
%
% DESCRIPTION  
% It defines the product combiner on a set of classifiers, e.g. 
% V=[V1,V2,V3] trained on the same classes, by selecting the class 
% which yields the highest value of the product of the classifier 
% outputs. This might also be used as A*[V1,V2,V3]*PRODC in which 
% A is a dataset to be classified.
% 
% If it is desired to operate on posterior probabilities, then the 
% input classifiers should be extended like V = V*CLASSC.
% 
% The base classifiers may be combined in a stacked way (operating
% in the same feature space) by V = [V1,V2,V3,...] or in a parallel
% way (operating in different feature spaces) by V = [V1;V2;V3;...].
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, VOTEC, MAXC, MINC, MEDIANC, MEANC, AVERAGEC, 
% STACKED, PARALLEL
%
% EXAMPLES
% See PREX_COMBINING.

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands 

% $Id: prodc.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function w = prodc(p1)

		type = 'prod';              % Define the operation processed by FIXEDCC.
	name = 'Product combiner';  % Define the name of the combiner.

	% This is a general procedure for all possible 
	% calls of fixed combiners handled by FIXEDCC
	if nargin == 0
		w = prmapping('fixedcc','combiner',{[],type,name});
	else
		w = fixedcc(p1,[],type,name);
	end

	if isa(w,'prmapping'),
		w = setname(w,name);
	end
	
	return
