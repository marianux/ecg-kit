%FEATSEL Selection of known features
%
%   [W,V] = FEATSEL(K,J)
%
% INPUT
%   K    Input dimensionality
%   J    Index vector of features to be selected
%
% OUTPUT
%   W    Mapping performing the feature selection
%   V    Mapping selecting the complementing features 
%
% DESCRIPTION
% This is a simple support routine that writes feature selection
% in terms of a mapping. If A is a K-dimensional dataset and J are
% the feature indices to be selected, then B = A*W does the same as
% B = A(:,J).
%
% The use of this routine is a mapping V computed for a lower dimensional
% subspace defined by J can now be defined by W = FEATSEL(K,J)*V as a 
% mapping in the original K-dimensional space.
%
% The selected features can be retrieved by W.DATA or by +W.
% See below for various methods to perform feature selection.
%
% SEE ALSO
% MAPPINGS, DATASETS, FEATEVAL, FEATSELF, FEATSELLR,
% FEATSELO, FEATSELB, FEATSELI, FEATSELP, FEATSELM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [w,v] = featsel(k,j)

if isa(k,'prdataset') | (isdouble(k) & (numel(k)>1))
%if isa(k,'prdataset') & ismapping(j)
	nodatafile(k);
	w = k(:,j);
  if nargout > 1
    n = [1:size(k,2)];
    n(j) = [];
    v = k(:,n);
  end
elseif isa(k,'prmapping')
	w = k(:,j);
  if nargout > 1
    n = [1:size(k,2)];
    n(j) = [];
    v = k(:,n);
  end
else
	if (any(j) > k | any(j) < 1)
		error('Features to be selected are not in proper range')
	end
	% w = prmapping(mfilename,'trained',j(:)',[],k,length(j));
	% There seems to be no need to make this mapping 'trained'.
	% A fixed mapping could be used more easily.
	% For historical consistency it is left like this (RD)
	% On second sight a fixed mapping is needed for handling datafiles
	% w = prmapping(mfilename,'trained',j(:)',[],k,length(j));
	w = prmapping(mfilename,'fixed',j(:)',[],k,length(j));
	% w = prmapping(mfilename,'combiner',j(:)',[],k,length(j));
	w = setname(w,'Feature Selection');
end
