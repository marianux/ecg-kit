%IM_MINF Minimum filter of images stored in a dataset (DIP_Image)
%
%	B = IM_MINF(A,SIZE,SHAPE)
%	B = A*IM_MINF([],SIZE,SHAPE)
%
% INPUT
%   A        Dataset with object images dataset (possibly multi-band)
%   SIZE     Filter width in pixels, default SIZE = 7
%   SHAPE    String with shape:'rectangular', 'elliptic', 'diamond'
%            Default: elliptic
%
% OUTPUT
%   B        Dataset with filtered images
%
% SEE ALSO
% DATASETS, DATAFILES, DIP_IMAGE, MINF

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_minf(a,size,shape)

		
	if nargin < 3 | isempty(shape), shape = 'elliptic'; end
	if nargin < 2 | isempty(size), size = 7; end
	
  if nargin < 1 | isempty(a)
    b = prmapping(mfilename,'fixed',{size,shape});
    b = setname(b,'Minimum filter');
	elseif isa(a,'prdataset') % allows datafiles too
		isobjim(a);
    b = filtim(a,mfilename,{size,shape});
  elseif isa(a,'double') | isa(a,'dip_image') % here we have a single image
		a = 1.0*dip_image(a);
		b = minf(a,size,shape);
	end
	
return
