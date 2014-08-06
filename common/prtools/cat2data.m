%CAT2DATA Create categorical dataset
%
%   A = CAT2DATA(DATA,LABELS)
%
% INPUT
%   DATA         Cell array or integer array with categorical data
%   LABELS       Desired labels, default: unlabeled
%
% OUTPUT
%   A            Dataset
%
% DESCRIPTION
% This routine creates a dataset of categorical data, either given as
% integers, or as a set of string arrays. If given the dataset is labeled
% by LABELS. In case DATA is a cell array of size [1,K], it is assumed that
% every cell contains a string array of M named categories, in which M is
% the number of objects. In case DATA is a cell array of size [M,1], it is
% assumed that every cell contains a string array of K named categories,
% one for every feature. Alternatively DATA can be a [M,K] cell array
% containing a single string per cell.
%
% The routine FEAT2LAB may be used to convert a categorical feature into a
% class label.
%
% A categorical dataset A can be combined with a dataset B representing a
% set of numeric features for the same objects by [A B].
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, PRDATASET, SETFEATDOM, GETFEATDOM, FEAT2LAB

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function a = cat2data(data,labels)

if nargin < 2, labels = []; end

if isdouble(data)
  [m,k] = size(data)
  a = prdataset(data,labels);
  featdom = cell(1,k);
  for j=1:k
    featdom{j} = unique(data(:,j))';
  end
elseif ischar(data)
  [ndat,lablist]=renumlab(data);
  a = prdataset(ndat,labels);
  featdom = {lablist};
elseif iscell(data) && size(data,1) == 1
  % data is given feature wise: as char arrays over columns
  k = numel(data);
  m = size(data{1},1);
  ndat = zeros(m,k);
  featdom = cell(1,k);
  for j=1:k
    [ndatj,lablistj]=renumlab(data{j});
    if size(ndat,1) ~= m
      error('Each set of nominal features should be based on the same number of objects')
    end
    ndat(:,j) = ndatj;
    featdom{j} = lablistj;
  end
  a = prdataset(ndat,labels);
elseif iscell(data) && size(data,2) == 1
  % data is given object wise: as char arrays over rows, Let's transpose it
  m = numel(data);
  k = size(data{1},1);
  data2 = cell(m,k);
  data1 = cell(1,k);
  for i=1:m
    data2(i,:) = cellstr(data{i});
  end
  for j=1:k
    data1{j} = char(data2(:,j));
  end
  a = feval(mfilename,data1,labels);
  return
elseif iscell(data)
  % data is given as a 2D cell array with a single string per cell
  [m,k] = size(data);
  data1 = cell(1,k);
  for j=1:k
    data1{j} = char(data(:,j));
  end
  a = feval(mfilename,data1,labels);
  return
else
  error('Data should be given as doubles or as a char cell array')
end
a = setfeatdom(a,featdom);

return