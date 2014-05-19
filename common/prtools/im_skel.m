%IM_SKEL Skeleton of binary images stored in a dataset (DIP_Image)
%
%	B = IM_SKEL(A)
%	B = A*IM_SKEL
%
% INPUT
%   A        Dataset with binary object images dataset 
%
% OUTPUT
%   B        Dataset with skeleton images
%
% SEE ALSO
% DATASETS, DATAFILES, DIP_IMAGE, BSKELETON

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_berosion(a)

		
  if nargin < 1 | isempty(a)
    b = prmapping(mfilename,'fixed');
    b = setname(b,'Image skeleton');
	elseif isa(a,'prdataset') % allows datafiles too
		isobjim(a);
    b = filtim(a,mfilename);
  elseif isa(a,'double') | isa(a,'dip_image') % here we have a single image
		a = dip_image(a,'bin');
		b = bskeleton(a,0,'natural');
	end
	
return
