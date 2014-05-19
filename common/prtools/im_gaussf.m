%IM_GAUSSF Gaussian filter of images stored in a dataset (DIPImage)
%
%	B = IM_GAUSSF(A,S)
%	B = A*IM_GAUSSF([],S)
%
% INPUT
%   A        Dataset with object images dataset (possibly multi-band)
%   S        Desired standard deviation for filter, default S = 1
%
% OUTPUT
%   B        Dataset with Gaussian filtered images
%
% DESCRIPTION
% All, possibly multi-band, 2D images in A are Gaussian filtered using the
% DIPImage command GAUSSF.
% In case DIPImage is not available, IM_GAUSS may be used.
%
% SEE ALSO
% DATASETS, DATAFILES, DIPIMAGE, GAUSSF, IM_GAUSS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_gaussf(a,s)

		
	if ~isdipimage
		error('DipImage not available, use im_gauss instead of im_gaussf')
	end
	
	if nargin < 2 | isempty(s), s = 1; end
	
  if nargin < 1 | isempty(a)
    b = prmapping(mfilename,'fixed',s);
    b = setname(b,'Gaussian filter');
	elseif isa(a,'prdataset') % allows datafiles too
		isobjim(a);
    b = filtim(a,mfilename,s);
  elseif isa(a,'double') | isa(a,'dip_image') % here we have a single image
		n = size(a,3);
		b = zeros(size(a));
		for j=1:n
			aa = 1.0*dip_image(a(:,:,j));
    	b(:,:,j) = double(gaussf(aa,s));
		end
	end
	
return
