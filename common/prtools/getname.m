%GETNAME High level routine for finding names of datasets and classifiers
% 
%   NAME = getname(A,N)
% 
% INPUT
%   A    Dataset, mapping or cell array of datasets or mappings
%   N    Number of characters in NAME (default: all)
% 
% OUTPUT
%   NAME  Dataset name or cell array of names
% 
% DESCRIPTION
% If N given, the return string has exactly N characters. This is done by
% truncation or by padding with blanks. This is useful for display purposes.
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function name = getname(c,n)

  if nargin < 2, n = []; end
  if iscell(c)
    name = cell(1,numel(c));
    for j=1:numel(c)
      name{j} = getname(c{j},n);
    end
  elseif isa(c,'prdataset') || ismapping(c)
    % this will go to dataset or mapping getname
    name = getname(c,n);
  else
    error('Illegal data type')
  end

return