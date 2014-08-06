%COL2GRAY Fixed mapping converting multi-band images to single band images
%
% COL2GRAY is identical to IM_GRAY, preserved for historical consistency.

function b = col2gray(varargin)
	 
  b = im_gray(varargin{:});
  
return