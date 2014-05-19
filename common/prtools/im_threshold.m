%IM_THRESHOLD Threshold images stored in a dataset (DIP_Image)
%
%	B = IM_THRESHOLD(A,TYPE,PAR,INV)
%	B = A*IM_THRESHOLD([],TYPE,PAR,INV)
%
% INPUT
%   A        Dataset with object images (possibly multi-band)
%   TYPE     Type of procedure, see below
%   PAR      Related parameter
%   INV      If INV = 1, result inverted, default INV = 0.
%
% OUTPUT
%   B        Dataset with thresholded images
%
% DESCRIPTION
% The following procedures are supported (TYPE)
% 'isodata':    Thresholding using the Isodata algorithm,
%               for more options (mask image, several thresholds)
%               see dip_isodatathreshold. (default)
% 'triangle':   Thresholding using chord method
%               (a.k.a. skewed bi-modality, maximum distance to triangle)
%               by Zack, Rogers and Latt (1973).
% 'background': Thresholding using unimodal background-symmetry method.
% 'fixed':      Thresholding at a fixed value.
% 'double':     Thresholding between two fixed values.
% 'volume':     Thresholding to obtain a given volume fraction.
% 'hysteresis': From the binary image (in>low) only those regions are
%               selected for which at least one pixel is (in>high)
%
% The following parameters are related to these procedures (PAR):
% ('background'): Distance to the peak where we cut-off, in
%               terms of the half-width at half the maximum.
%               Inf selects the default value, which is 2.
% ('fixed'):    Threshold value. Inf means halfway between
%               minimum and maximum value.
% ('double'):   Two threshold values. Inf means min+[1/3,2/3]*(max-min).
% ('volume'):   Parameter = the volume fraction (Inf means 0.5)
% ('hysteresis'):Two values: low, high (see above)
%
% SEE ALSO
% DATASETS, DATAFILES, DIP_IMAGE, THRESHOLD

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_threshold(a,type,par,inv)

		
	if nargin < 4 | isempty(inv), inv = 0; end
	if nargin < 3 | isempty(par), par = inf; end
	if nargin < 2 | isempty(type), type = 'isodata'; end
	
  if nargin < 1 | isempty(a)
    b = prmapping(mfilename,'fixed',{type,par,inv});
    b = setname(b,'Image threshold');
	elseif isa(a,'prdataset') % allows datafiles too
		isobjim(a);
    b = filtim(a,mfilename,{type,par,inv});
		b = setfeatsize(b,getfeatsize(a));
  elseif isa(a,'double') | isa(a,'dip_image') % here we have a single image
		a = 1.0*dip_image(a);
		b = threshold(a,type,par);
		if inv
			b = 1-b;
		end
	end

return