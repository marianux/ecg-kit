%IM_FFT Fixed mapping for 2D FFT
%
%	F = IM_FFT(A,POST1,POST2,...)
% F = A*IM_FFT([],POST1,POST2,...)
% F = A*IM_FFT(POST1,POST2,...)
%
% INPUT
%   A        Dataset with object images (possibly multi-band)
%
% OUTPUT
%   F        Dataset with FFT images
%   POST1    Characterstring defining postprocessing by FILTIM
%   POST2    Characterstring defining more postprocessing by FILTIM
%
% DESCRIPTION
% This routine applies the 2D FFT2 routine to all images in A. The result
% is complex. See the below example for post processing by FILTIM. It shows
% how it can be integrated in the call to IM_FFT.
%
% EXAMPLE
% delfigs
% prdatafiles;        % make sure PRDATAFILES is in the path
% a = delft_idb;      % load the Delft Image DataBase
% x = a*gendat(4);    % load a subset of 4 images
% % show the centralized log of the powerspectra of the gray images
% show(x*im_gray*im_fft*filtim('abs')*filtim('fftshift')*filtim('log'));
% % post-processing can be included in im_fft
% figure; show(x*im_gray*im_fft('abs','fftshift','log'));
% figure; show(x);    % show the originals
% showfigs
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES, FFT2, FILTIM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function b = im_fft(varargin)

  if nargin==0, varargin = {[]}; end
	varargin = shiftargin(varargin,'char');
  a = varargin{1};
	if isempty(a)
		b = prmapping(mfilename,'fixed',varargin(2:end));
		b = setname(b,'Image FFT');
	elseif isa(a,'prdataset') % allows datafiles too
		isobjim(a);
		b = filtim(a,mfilename,varargin(2:end));
		b = setfeatsize(b,getfeatsize(a));
	elseif isa(a,'double') || isa(a,'dip_image') % here we have a single image
		a = double(a);
		b = fft2(a);
		if numel(varargin) > 1
			for j=2:numel(varargin)
				b = filtim(b,varargin{j});
			end
    end
	end

return