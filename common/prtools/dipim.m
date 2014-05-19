%DIPIM Run any DIPimage command with one input image (DIP_Image)
%
%	B = DIPIM(A,COMMAND,PAR1,PAR2, ...)
%	B = A*DIPIM([],COMMAND,PAR1,PAR2, ...)
%
% INPUT
%   A        Dataset or datafiles with object images
%   COMMAND  String with name of a DIPimage command
%   PAR1     Additional parameters needed for COMMAND
%
% OUTPUT
%   B        Datasetor datafile with results
%
% DESCRIPTION
% This is a general routine to run DIPimage commands on the images stored
% in a dataset or datafile. It is thereby appropriate for grey value 
% operations. For logical and morphological operations DIPBIN should be
% used. See DIPIMAGE for a list of all possible commands. A few suggections:
%     stretch              - Grey-value stretching
%     hist_equalize        - Histogram equalization
%     lut                  - Look-up Table (with interpolation)
%     convolve             - General convolution filter
%     gaussf               - Gaussian filter
%     unif                 - Uniform filter
%     maxf                 - Maximum filter
%     minf                 - Minimum filter
%     medif                - Median filter
%     percf                - Percentile filter
%     varif                - Variance filter
%     gabor                - Gabor filter
%     derivative           - Derivative filters
%     dx                   - First derivative in the X-direction
%     dy                   - First derivative in the Y-direction
%     dz                   - First derivative in the Z-direction
%     laplace              - Laplace operator
%     laplace_plus_dgg     - Laplace + Dgg
%     laplace_min_dgg      - Laplace - Dgg
%     hessian              - Hessian matrix of an image
%     prewittf             - Prewitt derivative filter
%     sobelf               - Sobel derivative filter
%     threshold            - Thresholding
%     hist2image           - Backmaps a 2D histogram ROI to the images
%     minima               - Detect local minima
%     maxima               - Detect local maxima
%     watershed            - Watershed
% 
% SEE ALSO
% DATASETS, DATAFILES, DIP_IMAGE

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = dipim(a,command,varargin)

		
	if nargin < 2
    error('No DIPimage command supplied')
  end
	
  if isempty(a)
    b = prmapping(mfilename,'fixed',{command varargin{:}});
    b = setname(b,[command '(DIPimage)']);
	elseif isa(a,'prdataset') % allows datafiles too
		isobjim(a);
    b = filtim(a,mfilename,{command varargin{:}});
  elseif isa(a,'double') % here we have a single image
		a = dip_image(a);
    b = runcommand(command,a,varargin{:});
  elseif isa(a,'dip_image') % here we have a single DIP image
    b = runcommand(command,a,varargin{:});
  else
    error('Illegal call')
  end
  
 return
  
  function b = runcommand(command,a,varargin)
    
  try
    b = feval(command,a,varargin{:});
  catch ME
    if strcmp(ME.message, 'Argument # 1: image data type not supported')
      error('Wrong routine for this DIPimage coomand. Try DIPBIN.')
    end
    rethrow(ME);
  end
    
return
