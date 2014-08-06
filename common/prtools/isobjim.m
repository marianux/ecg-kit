%ISOBJIM test if the dataset contains objects that are images
%
%  N = ISOBJIM(A)
%      ISOBJIM(A)
%
% INPUT
%  A  input dataset
%
% OUTPUT
%  N  logical value
%
% DESCRIPTION
% True if the dataset contains objects that are images. If no output is required,
% false outputs are turned into errors. This may be used for assertion.
%
% In case the objects in A should be interpreted as 1-dimensional object
% images, use  A = SETFEATSIZE(A,[1 GETFEATSIZE(A)]) for conversion.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% ISDATASET, ISDATAIM, SETFEATSIZE, GETFEATSIZE

% $Id: isobjim.m,v 1.5 2009/01/31 15:43:10 duin Exp $

function n = isobjim(a)

		
	n = (isa(a,'prdataset') & length(a.featsize) > 1) | isdatafile(a);

	% generate error if input is not a dataset with image data with
   % pixels being features AND no output is requested (assertion)

	if nargout == 0 & n == 0
		error([newline '---- Dataset with object images expected -----'])
	end

return
