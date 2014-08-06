%IM_INVERT Fixed mapping for images inversion
%
% A = IM_INVERT(A)
% A = A*IM_INVERT
%
% DESCRIPTION
% Inverts image A by subtracting it from its maximum. Note that binary
% images can better be inverted by A = ~A. In that case also A = 1-A can
% be done.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES

% Copyright: D. de Ridder, R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_invert(a)
		
  if nargin < 1 | isempty(a)
    b = prmapping(mfilename,'fixed');
    b = setname(b,'Image inverse');
  elseif isa(a,'prdataset') % allows datafiles too
    isobjim(a);
    b = filtim(a,mfilename);
  elseif isa(a,'double') | isa(a,'dip_image') % here we have a single image
		if isa(a,'dip_image'), a = double(a); end
		b = max(max(max(a)))-a;
  else
    error('Illegal input')
	end
	
return

