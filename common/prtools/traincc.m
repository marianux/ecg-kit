%TRAINCC Train combining classifier if needed
% 
%   W = TRAINCC(A,W,CCLASSF)
% 
% INPUT
%   A        Training dataset
%   W        A set of classifiers to be combined  
%   CCLASSF  Combining classifier
%
% OUTPUT
%   B        Combined classifier mapping
%
% DESCRIPTION  
% The combining classifier CCLASSF is trained by the dataset A*W, if
% training is needed. W is typically a set of stacked (operating in the same
% feature space) or parallel (operating in different feature spaces;
% performed one after another) classifiers to be combined. E.g. if V1, V2
% and V3 are base classifiers, then V = [V1,V2,V3,...] is a stacked
% classifier and V = [V1;V2;V3;...] is a parallel one. If CCLASSF is one of
% the fixed combining rules like MAXC, then training is skipped.
%
% This routine is typically called by combining classifier schemes like
% BAGGINGC and BOOSTINGC.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, STACKED, PARALLEL, BAGGINGC

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: traincc.m,v 1.4 2010/06/25 07:55:34 duin Exp $

function w = traincc(a,w,cclassf)

		if (~ismapping(cclassf))
		error('Combining classifier is an unknown mapping.')
	end

	% If CCLASSF is already a combining classifier, just apply it. Otherwise,
	% train it using A*W.
	if isuntrained(w)
		w = a*w;  % train base classifiers
	end

	w = w*cclassf;
	if isuntrained(w)
		w = a*w;  % train combiner when needed
	end

	w = setcost(w,a);
	
return
