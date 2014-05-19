%DIPIMAGECHECK  Check whether the DIP_IMAGE toolbox is in the path

function dipimagecheck
	
	if exist('dip_image','file') ~= 3
		error([newline 'The DIPIMAGE package is needed and not found.' ...
      newline 'Please download and install it from http://www.diplib.org/'])
	end
		