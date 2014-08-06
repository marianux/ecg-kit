%DATFILT Filtering of dataset images (deprecated)
%
%   B = DATFILT(A,F)
%
% INPUT
%   A  Dataset with image data
%   F  Matrix with the convolution mask
%
% OUTPUT
%   B  Dataset containing all the images after filtering
%
% DESCRIPTION
% All images stored in the dataset A are horizontally and vertically
% convoluted by the 1-dimensional filter F. A uniform N*N filter is, 
% thereby, realized by DATFILT(A,ONES(1,N)/N).
% 
% This command is deprecated, use IMFILT instead.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, IM2OBJ, DATA2IM, IM2FEAT, DATGAUSS, DATAIM

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: datfilt.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function a = datfilt(a,f)
		
	[m,k] = getsize(a);
	n  = length(f);
	nn = floor(n/2);
	im = data2im(a);
	[imheight,imwidth,nim] = size(im);

	for i=1:nim

		% Add a border with NN pixels, set the border to
		% the mirrored original values (private function).
		c = bord(im(:,:,i),NaN,nn);
    c = conv2(f,f,c,'same');
		im(:,:,i) = resize(c,nn,imheight,imwidth);
	end

	if (isfeatim(a))
		a = setdata(a,im2feat(im),getfeatlab(a));
	else
		a = setdata(a,im2obj(im),getfeatlab(a));
	end
	
return;
