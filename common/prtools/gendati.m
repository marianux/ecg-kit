%GENDATI Create dataset from randomly selected windows in a given image
%
%  A = GENDATI(IMAGE,WSIZE,N,LABEL)
%
% INPUT
%   IMAGE  - Image of any dimensionality
%   WSIZE  - Vector with size of the window
%   N      - Number windows to be generated
%   LABEL  - Optional string or number with label for all objects
%
% OUTPUT
%   A      - Dataset of N objects and prod(WSIZE) features
%
% DESCRIPTION
% Windows of the specified size are arbitrarily positioned in the
% image, converted to a row vector and stored in the dataset A.
%
% If specified, all objects have label LABEL. Otherwise they are
% unlabeled.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function a = gendati(im,ws,m,label)

if nargin < 4, label = []; end
if nargin < 3, m = 100; end

n = length(ws);
imsize = size(im);
if length(imsize) ~= n
	error('Window size should have same number of components as image size')
end

if any(ws > imsize)
	error('Window size should not be larger than image size')
end

we = prod(ws);
n = length(ws);
window_offset = [0:ws(1)-1]';
for j=2:n
	wo = window_offset;
	for i=1:ws(j)-1
		wo = wo + prod(imsize(1:j-1));
		window_offset = cat(j,window_offset,wo);
	end
end
window_offset = window_offset(:)';

N = cell(1,n);
for j=1:n
	R = imsize(j)-ws(j)+1;
	N(j) = {ceil(rand(m,1)*R)};
end
N = sub2ind(imsize,N{:});

a = zeros(m,prod(ws));
for i=1:m
	a(i,:) = im(N(i)+window_offset);
end

a = prdataset(a);
a = setfeatsize(a,ws);
if ~isempty(label)
	a = setlabels(a,label);
end
