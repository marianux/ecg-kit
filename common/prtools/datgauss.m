%DATGAUSS Apply Gaussian filter on images in a dataset (deprecated)
%
%  B = DATGAUSS(A,SIGMA)
%
% INPUT
%  A        Dataset containing images
%  SIGMA    Standard deviation of Gaussian filter (default 1)
%
% OUTPUT
%  B       Dataset with filtered images
%
% DESCRIPTION
% All images stored as objects (rows) or as features (columns) of dataset A
% are filtered with a Gaussian filter with standard deviation SIGMA and 
% stored in dataset B. Image borders are mirrored before filtering.
%
% This command is deprecated, use IM_GAUSS instead.
%
% SEE ALSO
% DATASETS, DATAIM, IM2OBJ, IM2FEAT, DATUNIF, IM_GAUSS

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: datgauss.m,v 1.5 2009/01/04 21:04:44 duin Exp $

function a = datgauss(a,sigma)

		
	if nargin < 2 | isempty(sigma)
		sigma = 1;
	end
	if nargin < 1 | isempty(a)
		a = prmapping(mfilename,'fixed',{sigma});
		a = setname(a,'Gaussian filter');
		return
	end
	
	if ismapping(sigma)
		w = sigma;
		sigma = w.data{1};
	end

	% Convert dataset to image or image array.
 	im = data2im(a); [imheight,imwidth,nim] = size(im);
	
	% Check argument.
	sigma = sigma(:);
	if (length(sigma) ~= 1) & (length(sigma) ~= nim)
		error('incorrect mumber of standard deviations specified')
	end

	% Create filter(s): 1D Gaussian.
	bordersize = ceil(2*sigma); filtersize = 2*bordersize + 1;	
	filter = exp(-repmat((([1:filtersize] - bordersize - 1).^2),length(sigma),1) ...
								./repmat((2.*sigma.*sigma),1,filtersize));		% Gaussian(s).
	filter = filter ./ repmat(sum(filter,2),1,filtersize);			% Normalize.

	% Replicate filter if just a single SIGMA was specified.
	if (length(sigma) == 1)
		bordersize = repmat(bordersize,nim,1);
		filter     = repmat(filter,nim,1);
	end

	% Process all images...
	if exist('gaussf') == 2 
		for i=1:nim      
			im(:,:,i) = double(gaussf(1*im(:,:,i),sigma)); 
		end
	else
		for i = 1:nim
			out = bord(im(:,:,i),NaN,bordersize(i));           % Add mirrored border.
			out = conv2(filter(i,:),filter(i,:),out,'same');   % Convolve with filter.
			im(:,:,i) = resize(out,bordersize(i),imheight,imwidth);
																								% Crop back to original size.
		end
	end
	
  % Place filtered images back in dataset.
	im = squeeze(im);
	if (isfeatim(a))
		a = setdata(a,im2feat(im),getfeatlab(a));
	else
		a = setdata(a,im2obj(im,getfeatsize(a)),getfeatlab(a));
	end

return
