%IM_STRETCH Contrast stretching of images stored in a dataset (DIP_Image)
%
%	B = IM_STRETCH(A,LOW,HIGH,MIN,MAX)
%	B = A*IM_STRETCH([],LOW,HIGH,MIN,MAX)
%
% INPUT
%   A        Dataset with object images dataset (possibly multi-band)
%   LOW      Lower percentile (default 0)
%   HIGH     Highest percentile (default 100)
%   MIN      Miniumum (default 0)
%   MAX      Maximum (default 255)
%
% OUTPUT
%   B        Dataset with filtered images
%
% SEE ALSO
% DATASETS, DATAFILES, DIP_IMAGE, STRETCH

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_stretch(a,low,high,min,max)

		
	if nargin < 5 | isempty(max), max = 255; end
	if nargin < 4 | isempty(min), min = 0; end
	if nargin < 3 | isempty(high), high = 100; end
	if nargin < 2 | isempty(low), low = 0; end
	
  if nargin < 1 | isempty(a)
    b = prmapping(mfilename,'fixed',{low,high,min,max});
    b = setname(b,'Image stretch');
	elseif isa(a,'prdataset') % allows datafiles too
		isobjim(a);
    b = filtim(a,mfilename,{low,high,min,max});
  elseif isa(a,'double') | isa(a,'dip_image') % here we have a single image
		a = 1.0*dip_image(a);
		b = stretch(a,low,high,min,max);
	end
	
return
