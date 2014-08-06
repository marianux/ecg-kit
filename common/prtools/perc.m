%PERC Percentile combining classifier
% 
% 	W = PERC(V,P)
% 	W = V*PERC([],P)
% 
% INPUT
%   V     Set of classifiers
%   P     Percentile, 0 <= P <= 100
%
% OUTPUT
%   W    Percentile combining classifier on V
%
% DESCRIPTION
% If V = [V1,V2,V3, ... ] is a set of classifiers trained on the 
% same classes and W is the percentile combiner: it selects the class 
% defined by the percentile of the outputs of the input classifiers. This 
% might also be used as A*[V1,V2,V3]*PERC([],P) in which A is a dataset to 
% be classified. 
%
% PERC([],0)   is equal to MINC
% PERC([],50)  is equal to MEDIANC
% PERC([],100) is equal to MAXC
% 
% If it is desired to operate on posterior probabilities then the 
% input classifiers should be extended like V = V*CLASSC;
%
% The base classifiers may be combined in a stacked way (operating
% in the same feature space by V = [V1,V2,V3, ... ] or in a parallel
% way (operating in different feature spaces) by V = [V1;V2;V3; ... ]
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, VOTEC, MAXC, MINC, MEANC, MEDIANC, PRODC,
% AVERAGEC, STACKED, PARALLEL
%
% EXAMPLES
% See PREX_COMBINING

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function w = perc(p1,par)

  if nargin < 2 | par < 0 | par > 100
		error('Percentile between 0 and 100 should be defined for percentile combiner')
	end 

	type = 'perc'; % define the operation processed by FIXEDCC.

	% define the name of the combiner. 
	% this is the general procedure for all possible calls of fixed combiners
	% handled by FIXEDCC
	name = 'Percentile combiner'; 

	if nargin == 0
		w = prmapping('fixedcc','combiner',{[],type,name,par});
	else
		w = fixedcc(p1,[],type,name,par);
	end

	if isa(w,'prmapping')
		w = setname(w,name);
	end

return
	
