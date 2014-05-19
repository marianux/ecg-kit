%GETWINDOWS Get pixel feature vectors around given pixels in image dataset
%
%		L = GETWINDOWS(A,INDEX,WSIZE,INCLUDE)
%		L = GETWINDOWS(A,[ROW,COL],WSIZE,INCLUDE)
%
% INPUT
%   A       Dataset containing feature images
%   INDEX   Index vector of target pixels in the images (Objects in A)
%   ROW     Column vector of row-coordinates of target pixels
%   COL     Column vector of column-coordinates of target pixels
%   WSIZE   Desired size of rectangular window around target pixels
%   INCLUDE Flag (0/1), indicating whether target pixels should be included
%           (1,default), or not (0) in result.
% OUTPUT
%   L       Index in A of window pixels
%
% DESCRIPTION
% This routine generates all objects in a dataset constructed by image
% features that are in a window of size WSIZE around the target pixels
% given by INDEX, or by [ROW,COL]. If WSIZE is omitted or empty ([],
% default), just the 4 4-conntected neighbors of the target pixels are
% returned. L points to the window pixels, such that A(L,:) is a dataset
% of the corresponding objects.
%
% SEE ALSO
% DATASETS, IM2FEAT

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function L = getwindows(a,index,wsize,include) 

if nargin < 4, includeflag = 1; end
if nargin < 3, wsize = []; end

isdataset(a);
isfeatim(a);

imsize = getobjsize(a);

if size(index,2) == 1
	[row,col] = ind2sub(imsize,index);
elseif size(index,2) == 2
	row = index(:,1);
	col = index(:,2);
end
n = length(row);

imsize = getobjsize(a);

if ~isempty(wsize) % rectangular neighborhoods
	if length(wsize) == 1 % square windows
		row1 = ceil(wsize/2)-1; col1 = row1; % # pixels above / left of target pixel
		row2 = floor(wsize/2);  col2 = row2; % # pixels below / right of target pixel
	elseif length(wsize) == 2 % possibly non-square windows
		row1 = ceil(wsize(1)/2)-1;
		row2 = floor(wsize(1)/2);
		col1 = ceil(wsize(2)/2)-1;
		col2 = floor(wsize(2)/2);
	else
		error('Window size should be 1D or 2D')
	end
	[R,C] = meshgrid([-row1:row2],[-col1:col2]);
	k = length(R(:));
	R = repmat(row(:),1,k)+repmat(R(:)',length(row),1);
	C = repmat(col(:),1,k)+repmat(C(:)',length(col),1);
else
	R = repmat(row(:),1,5)+repmat([0 -1 0 1 0],length(row),1);
	C = repmat(col(:),1,5)+repmat([-1 0 0 0 1],length(col),1);
end
JR = [find(R <= 0); find(R > imsize(1))];
R(JR) = []; C(JR) = [];
JC = [find(C <= 0); find(C > imsize(2))];
R(JC) = []; C(JC) = [];
L = sub2ind(imsize,R,C);
%b(L) = ones(length(L),1);

if include % we are done
	L = unique(L);
else       % remove given objects
	b = zeros(imsize);             % create an image of the right size
	b(L) = ones(length(L),1);      % flag the objects we found.
	Z = sub2ind(imsize,row,col);   % for all original objects,
	b(Z) = zeros(length(Z),1);     % remove flags
	L = find(b > 0);               % and see what is left
end
%b = a(L,:);

L = unique(L);
J = [find(L<1); find(L>prod(imsize))];
L(J) = [];
%b = a(L,:);

