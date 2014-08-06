%ISFEATIM
%
% N = ISFEATIM(A);
%
% INPUT
%   A   Input dataset
%
% OUTPUT
%   N   1/0 if dataset A does/doesn't contain images
%
% DESCRIPTION
% True if dataset contains features that are images.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% ISDATASET, ISMAPPING, ISDATAIM

% $Id: isfeatim.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function n = isfeatim(a)
			% When the field objsize contains a vector instead of a scalar, the
	% features inside the dataset are images:
	n = isa(a,'prdataset') & length(a.objsize) > 1;

	if (nargout == 0) & (n == 0)
		error([newline '---- Dataset with feature images expected -----'])
	end

return;
