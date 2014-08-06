%FEATSEL Fixed mapping for selecting given features
%
%   [W,V] = FEATSEL(K,J)
%    B    = A*FEATSEL(J)
%    C    = A*FEATSEL(J,true)
%
% INPUT
%   K    Input dimensionality
%   J    Index vector of features to be selected
%   A    Dataset
%
% OUTPUT
%   W    Mapping performing the feature selection
%   V    Mapping selecting the complementing features
%   B    Dataset, A(J,:)
%   C    Dataset, A with features J removed
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
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, FEATEVAL, FEATSELF, FEATSELLR,
% FEATSELO, FEATSELB, FEATSELI, FEATSELP, FEATSELM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [w,v] = featsel(varargin)

if nargin == 2 && isscalar(varargin{1})
  % preserving old type of call: [w,v] = featsel(k,j)
  [k,j] = deal(varargin{:});
  R = [1:k];
  R(j) = [];
	w = prmapping(mfilename,'fixed',j(:)',[],0,length(j));
	v = prmapping(mfilename,'fixed',R(:)',[],0,length(R));
	w = setname(w,'Feature Selection');
	v = setname(v,'Feature Selection');
  return
end
  
% new type of call: b = a*featsel(j)
argin = shiftargin(varargin,'vector');
argin = setdefaults(argin,[],[],false);
if mapping_task(argin,'definition')
  w = define_mapping(argin,'fixed');
  w = setname(w,'Feature Selection');
else
  [a,j,inv] = deal(argin{:});
  if ismapping(j)
    j = +j;
  end
  if inv
    R = [1:size(a,2)];
    R(j) = [];
  else
    R = j;
  end
  w = a(:,R);
end
  
