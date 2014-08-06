%FEAT2LAB  Label dataset by one of its features and remove this feature
%
%   B = FEAT2LAB(A,N)
%
% INPUT
%   A   Dataset
%   N   Integer, pointing to feature to be used as label
%
% OUTPUT
%   B   Dataset, feature N is removed and used for labeling
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, SETFEATDOM, GETFEATDOM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function a = feat2lab(a,n)

isdataset(a);
k = size(a,2);

if (n<1) | (n>k)
  error('Desired feature not in range')
end

lablist = getfeatdom(a,n);
nlab = +a(:,n);
a(:,n) = [];

if ~isempty(lablist) & ischar(lablist{1})
  a = setlablist(a,lablist{1});
  a = setnlab(a,nlab);
else
  a = setlabels(a,nlab);
end