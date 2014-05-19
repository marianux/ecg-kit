%IM_FFT 2D FFT of all images in dataset
%
%	F = IM_FFT(A)
%
% INPUT
%   A        Dataset with object images (possibly multi-band)
%
% OUTPUT
%   F        Dataset with FFT images
%
% SEE ALSO
% DATASETS, DATAFILES

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_fft(a,varargin)

		
	if nargin < 1 | isempty(a)
		b = prmapping(mfilename,'fixed',varargin);
		b = setname(b,'Image FFT');
	elseif isa(a,'prdataset') % allows datafiles too
		isobjim(a);
		b = filtim(a,mfilename,varargin);
		b = setfeatsize(b,getfeatsize(a));
	elseif isa(a,'double') | isa(a,'dip_image') % here we have a single image
		a = double(a);
		b = fft2(a);
		if nargin > 1
			for j=1:nargin-1
				b = filtim(b,varargin{j});
			end
		end
			
	end

return