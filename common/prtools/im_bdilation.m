%IM_BDILATION Binary dilation of images stored in a dataset (DIP_Image)
%
%	B = IM_BDILATION(A,N,CONNECTIVITY,EDGE_CONDITION)
%	B = A*IM_BDILATION([],N,CONNECTIVITY,EDGE_CONDITION)
%
% INPUT
%   A        Dataset with binary object images dataset (possibly multi-band)
%   N        Number of iterations (default 1)
%   CONNECTIVITY    See BDILATION
%   EDGE_CONDITION  Value of edge, default 0
%
% OUTPUT
%   B        Dataset with dilated images
%
% SEE ALSO
% DATASETS, DATAFILES, DIP_IMAGE, BDILATION

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

function b = im_bdilation(a,n,connect,edgecon)

		
	if nargin < 4 | isempty(edgecon), edgecon = 0; end
	if nargin < 3 | isempty(connect), connect = -2; end
	if nargin < 2 | isempty(n), n = 1; end
	
  if nargin < 1 | isempty(a)
    b = prmapping(mfilename,'fixed',{n,connect,edgecon});
    b = setname(b,'Image dilation');
	elseif isa(a,'prdataset') % allows datafiles too
		isobjim(a);
    b = filtim(a,mfilename,{n,connect,edgecon});
  elseif isa(a,'double') | isa(a,'dip_image') % here we have a single image
		a = dip_image(a,'bin');
		b = bdilation(a,n,connect,edgecon);
	end
	
return
