%IM_MEAN Fixed mapping computating the centers of gravity of images
%
%   B = IM_MEAN(A)
%   B = A*IM_MEAN
%
% INPUT
%   A        Dataset with object images dataset (possibly multi-band)
%
% OUTPUT
%   B        Dataset with centers-of-gravity replacing images 
%            (possibly multi-band). The first component is always measured
%            in the horizontal direction (X).
%
% DESCRIPTION
% Computes the centers of gravity of all images stored in A.
% This center is computed, relative to the top-left corner of the image.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_mean(a)

		
  if nargin < 1 | isempty(a)
    b = prmapping(mfilename,'fixed');
    b = setname(b,'Image centers');
	elseif isa(a,'prdataset') % allows datafiles too
		isobjim(a);
    b = filtim(a,mfilename);
  elseif isa(a,'double') % here we have a single image
		if isa(a,'dip_image'), a = double(a); end
		[m,n] = size(a);
		g = sum(sum(a));
		JX = repmat([1:n],m,1);
		JY = repmat([1:m]',1,n);
		mx = sum(sum(a.*JX))/g;
		my = sum(sum(a.*JY))/g;
		b = [mx/n my/m];
	end
	
return
