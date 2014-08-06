%DATAIM Image operation on dataset images (deprecated)
%
%    B = DATAIM(A,'IMAGE_COMMAND',PAR1,PAR2,....)
%
% INPUT
%    A               Dataset containing images
%    IMAGE_COMMAND   Function name
%    PAR1, ...       Optional parameters to IMAGE_COMMAND
%
% OUTPUT
%    B               Dataset containing images processed by IMAGE_COMMAND
%
% DESCRIPTION
% For each image stored in A (as either feature or object), performs
%
%    IMAGE_OUT = IMAGE_COMMAND(IMAGE_IN,PAR1,PAR2,....)
%
% and stores the result in dataset B.
%
% This command is deprecated, use FILTIM instead.
%
% EXAMPLES
% B = DATAIM (A,'CONV2',[-1 0 1; -1 0 1; -1 0 1],'same');
% Performs a convolution with a horizontal gradient filter (see CONV2).
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, IM2OBJ, DATA2IM, IM2FEAT, FILTIM

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: dataim.m,v 1.4 2007/03/22 08:53:12 duin Exp $

function b = dataim (a,command,varargin)

		isdataim(a);      % Assert that A contains images.
	im = data2im(a);  % Convert A to array of images.

	[m,k] = size(a); [height,width,n_im] = size(im);

	% Perform command on first image to check whether image size stays equal
	% (for feature images; the number of objects can change without problems).

	first = feval(command,im(:,:,1),varargin{:});
        first = double(first); % for DipLib users
	[newheight,newwidth] = size(first);
	if (isfeatim(a)) & ((newheight ~= height) | (newwidth ~= width))
		error('image size is not allowed to change')
	end

	% Process the remaining images.

	out = zeros(newheight,newwidth,n_im); out(:,:,1) = first;
	for i = 2:n_im
		out(:,:,i) = double(feval(command,im(:,:,i),varargin{:}));
	end

	% Convert the image array back into a dataset.

	if (isfeatim(a))
		b = setdata(a,im2feat(out),getfeatlab(a));   %images are stored in columns
	else
		b = setdata(a,im2obj(out));                  %images are stored in rows									
		b = setfeatsize(b,size(first));
	end

return
