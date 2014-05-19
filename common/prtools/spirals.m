%SPIRALS 194 objects with 2 features in 2 classes
%
%	A = SPIRALS
%	A = SPIRALS(M,N)
%
% Load the dataset in A, select the objects and features according to the
% index vectors M and N. This is one of the Spiral dataset implementations.
%
% See also DATASETS, PRDATASETS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function a = spirals(M,N);
if nargin < 2, N = []; end
if nargin < 1, M = []; end
a = prdataset('spirals',M,N);
a = setname(a,'Spirals');
