%IM_GAUSSF Fixed mapping for Gaussian filtering images (DIPImage)
%
%	B = IM_GAUSSF(A,S)
%	B = A*IM_GAUSSF([],S)
%	B = A*IM_GAUSSF(S)
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
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES, DIPIMAGE, GAUSSF, IM_GAUSS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_gaussf(varargin)

	argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],1);
  if mapping_task(argin,'definition')
    b = define_mapping(argin,'fixed');
    b = setname(b,'Gaussian filter');
  else
	  checktoolbox('diplib');
    [a,s] = deal(argin{:});	
    if isa(a,'prdataset') % allows datafiles too
      isobjim(a);
      b = filtim(a,mfilename,s);
    elseif isa(a,'double') | isa(a,'dip_image') % here we have a single image
      n = size(a,3);
      b = zeros(size(a));
      if checktoolbox('dipimage')
        for j=1:n
          aa = 1.0*dip_image(a(:,:,j));
          b(:,:,j) = double(gaussf(aa,s));
        end
      else
        diplibwarn
        b = im_gauss(a,s,s);
      end
    end
  end
	
return
