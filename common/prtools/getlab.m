%GETLAB Get labels of dataset or mapping
%
%  LABELS = GETLAB(A)
%  LABELS = GETLAB(W)
%
% INPUT
%  A      Dataset
%  W      Mapping
%
% OUTPUT
%  LABELS Labels
%
% DESCRIPTION
% Returns the labels of the objects in the dataset A or the feature labels
% assigned by the mapping W.
%
% If A (or W) is neither a dataset nor a mapping, a set of dummy
% labels is returned, one for each row in A. All these labels have the
% value NaN.
%
% Note that for datasets and mappings GETLAB() is effectively the same
% as GETLABELS().
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, GETLABELS

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: getlab.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function labels = getlab(a)

		if isa(a,'prdataset')
		labels = getlabels(a,'crisp');
	elseif isa(a,'prdatafile')
		labels = getlabels(a,'crisp');
	elseif isa(a,'prmapping')
		labels = getlabels(a);
	else
		labels = repmat(NaN,size(a,1),1);
	end

	return
