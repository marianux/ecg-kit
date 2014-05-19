%OBJ2FEAT Transform object images to feature images in dataset
%
%   B = OBJ2FEAT(A)
%
% INPUT
%   A     Dataset with object images, possible with multiple bands.
%
% OUTPUT
%   B     Dataset with features images.
%
% SEE ALSO
% DATASETS, IM2OBJ, IM2FEAT, DATA2IM, FEAT2OBJ

function b = obj2feat(a)

	
 	isdataset(a)
  isobjim(a);
  im = data2im(a);
  b = im2feat(im);
    
  
