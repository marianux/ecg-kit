%IM_LABEL Labeling of binary images stored in a dataset (DIP_Image)
%
%	B = IM_LABEL(A,CONNECTIVITY,MIN_SIZE,MAX_SIZE)
%	B = A*IM_LABEL([],CONNECTIVITY,MIN_SIZE,MAX_SIZE)
%
% INPUT
%   A        Dataset with binary object images dataset (possibly multi-band)
%   N        Number of iterations (default 1)
%   CONNECTIVITY    See LABEL
%   MIN_SIZE Minimum size of objects to be labeled (default 0)
%   MAX_SIZE Maximum size of objects to be labeled (default 0: all)
%
% OUTPUT
%   B        Dataset with labeled images
%
% SEE ALSO
% DATASETS, DATAFILES, DIP_IMAGE, LABEL

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_label(a,connect,minsize,maxsize)

		
	if nargin < 4 | isempty(maxsize), maxsize = 0; end
	if nargin < 3 | isempty(minsize), minsize = 0; end
	if nargin < 2 | isempty(connect), connect = 2; end
	
  if nargin < 1 | isempty(a)
    b = prmapping(mfilename,'fixed',{connect,minsize,maxsize});
    b = setname(b,'Image labeling');
	elseif isa(a,'prdataset') % allows datafiles too
		isobjim(a);
    b = filtim(a,mfilename,{connect,minsize,maxsize});
  elseif isa(a,'double') | isa(a,'dip_image') % here we have a single image
		a = dip_image(a,'bin');
		b = double(label(a,connect,minsize,maxsize));
	end
	
return
