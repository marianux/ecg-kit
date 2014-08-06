%CONCATM Fixed-cell mapping concatenating cell array of datasets or mappings
%
%   W = V*CONCATM
%   B = A*CONCATM
%
% INPUT
%   V  Cell array of mappings
%   A  Cell array of datasets
%
% OUTPUT
%   W  Concatenation of mappings: W = [V{:}]
%   A  Concatenation of datasets: B = [A{:}]
%
% DESCRIPTION
% This routine is an exception of the general treatment of cell arrays and
% mappings in PRTools. As a rule the elements of a cell array are processed
% one by one by a mapping. This routine takes the entire cell array and
% concatenates its elements.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function out = concatm(argin)

if nargin == 0
  out = prmapping(mfilename,'fixed_cell');
else
  if ~iscell(argin)
    error('Cell array expected')
  else
    out = [];
    for j=1:size(argin,1)
      out = [out;[argin{j,:}]];
    end
  end
end