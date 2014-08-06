%ISDATAIM Returns true if a dataset contains image objects or image features 
%
% 	N = ISDATAIM(A)
%
% INPUT
%		A 	Dataset
%
% OUTPUT
%		N	  Scalar: 1 if A contains images as objects or features, otherwise 0
%
% DESCRIPTION
% If no output argument is given, the function will produce an error if A does
% not contain image objects or features (i.e. it will act as an assertion).
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% ISOBJIM, ISFEATIM

% $Id: isdataim.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function n = isdataim (a)

		if (nargin < 1)
		error('Insufficient number of arguments.'); 
	end

	% The FEATSIZE and OBJSIZE fields of a dataset indicate whether it contains
	% images.
	
	n = (isa(a,'prdataset')) & ((length(a.featsize) > 1) | (length(a.objsize) > 1)) ;

	% Generate and error if the input is not a dataset with image data and
	% no output is requested (assertion).

	if (nargout == 0) & (n == 0)
		error('Dataset with images expected.')
	end

return
