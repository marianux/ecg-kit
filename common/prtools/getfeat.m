%GETFEAT Get feature labels of a dataset or a mapping
%
%	  LABELS = GETFEAT(A)
%	  LABELS = GETFEAT(W)
%
% INPUT
%   A,W    Dataset or mapping
%
% OUTPUT
%   LABELS Label vector with feature labels
%
% DESCRIPTION
% Returns the labels of the features in the dataset A or the labels
% assigned by the mapping W.
%
% Note that for a mapping W, getfeat(W) is effectively the same as GETLAB(W).
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>) 
% DATASETS, MAPPINGS, GETLAB

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: getfeat.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function labels = getfeat(par)

	if isdataset(par) | isdatafile(par)
		labels = getfeatlab(par);
	elseif isa(par,'prmapping')
		labels = getlabels(par);
	else
		error([newline 'Dataset or mapping expected.'])
	end
return;
