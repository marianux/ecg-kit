%ISSTACKED Test on stacked mapping
%
%   N = ISSTACKED(W)
%   ISSTACKED(W)
%
% INPUT
%   W  Mapping
%
% OUTPUT
%   N  Scalar, 1 if W is a stacked mapping, 0 otherwise
%
% DESCRIPTION
% Returns 1 for stacked mappings. If no output is requested, false outputs
% are turned into errors. This may be used for assertion.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% ISMAPPING, ISPARALLEL

% $Id: isstacked.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function n = isstacked(w)

		if (isa(w,'prmapping')) & (strcmp(w.mapping_file,'stacked'))
%			strcmp(w.mapping_file,'dyadicm')
		n = 1;
	else
		n = 0;
	end

	% Generate an error if the input is not a stacked mapping and no output 
	% is requested (assertion).

	if (nargout == 0) & (n == 0)
		error([newline '---- Stacked mapping expected -----'])
	end

return
