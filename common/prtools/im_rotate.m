%IM_ROTATE Rotate all images in dataset 
%
%	B = IM_ROTATE(A,ALF)
%
% INPUT
%   A        Dataset with object images (possibly multi-band)
%   ALF      Rotation angle (in radians), 
%            default: rotation to main axis
%
% OUTPUT
%   B        Dataset with rotated object images 
%
% SEE ALSO
% DATASETS, DATAFILES, DIP_IMAGE

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_rotate(a,alf)

		if nargin < 2, alf = []; end
	
  if nargin < 1 | isempty(a)
    b = prmapping(mfilename,'fixed',{alf});
    b = setname(b,'Image rotate');
	elseif isdataset(a)
		error('Command cannot be used for datasets as it may change image size')
	elseif isdatafile(a)
		isobjim(a);
    b = filtim(a,mfilename,{alf});
		b = setfeatsize(b,getfeatsize(a));
  elseif isa(a,'double') | isa(a,'dip_image') % here we have a single image

		a = double(a);

		if isempty(alf)
			m = im_moments(a,'central',[1 2 0; 1 0 2]');
			C = [m(2) m(1); m(1) m(3)];
			[E,D] = preig(C); [DD,ind] = sort(diag(D));
			alf = atan2(E(2,ind(1)),E(1,ind(1)))+pi/2;
		end

 		b = imrotate(a,alf*360/(2*pi),'bilinear','crop');
	
	end

return
