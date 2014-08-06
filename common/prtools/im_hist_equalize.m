%IM_HIST_EQUALIZE Fixed mapping for histogram equalization (DIP_Image)
%
%	B = IM_HIST_EQUALIZE(A)
%	B = A*IM_HIST_EQUALIZE
%
% INPUT
%   A        Dataset with object images dataset (possibly multi-band)
%
% OUTPUT
%   B        Dataset with filtered images
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES, DIP_IMAGE, HIST_EQUALIZE

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_hist_equalize(a)

  if nargin < 1 | isempty(a)
    b = prmapping(mfilename,'fixed');
    b = setname(b,'Image hist_equalize');
	elseif isa(a,'prdataset') % allows datafiles too
		isobjim(a);
    b = filtim(a,mfilename);
  elseif isa(a,'double') | isa(a,'dip_image') % here we have a single image
    if checktoolbox('dipimage')
      a = 1.0*dip_image(a);
      b = hist_equalize(a);
    else
      diplibwarn
      b = histeq(a);
    end
	end
	
return
