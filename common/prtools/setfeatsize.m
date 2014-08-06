%SETFEATSIZE Fixed mapping, (re)set the feature size of a dataset
%
%   A = SETFEATSIZE(A,FEATSIZE)
%   A = A*SETFEATSIZE([],FEATSIZE)
%   A = A*SETFEATSIZE(FEATSIZE)
%
% INPUT 
%   A         Dataset
%   FEATSIZE  Feature size vector, default: scalar with actual featsize
% 
% OUTPUT
%   A         Dataset
%
% DESCRIPTION
% By default the feature size of a dataset is its number of features, i.e.
% the number of columns in the DATA field of A. If the features are samples
% of a multi-dimensional data item, e.g. the pixels of an image, the
% original size of this data item may be stored in FEATSIZE. The product of
% all elements in FEATSIZE has to be equal to the number of columns in the
% DATA field of A.
%
% This routine is particulary useful as A = A*SETFEATSIZE to remove image
% size settings for the feature size.
%
% There is a low level version of this routine with the same name in the
% PRDATASET directory.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function out = setfeatsize(varargin)

  varargin = shiftargin(varargin,'double');
  argin = setdefaults(varargin,[],[]);
  if mapping_task(argin,'definition')
    out = define_mapping(argin,'fixed')
  else
    % this will not happen. prmap directly jumps to @prdataset/setfeatsize
    [a,featsize] = deal(argin{:});
    if isempty(featsize)
      featsize = prod(getfeatsize(a));
    end
    out = setfeatsize(a,featsize);
  end