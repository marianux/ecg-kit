%ISMAPPING Test whether the argument is a mapping
%
% 	N = ISMAPPING(W);
%
% INPUT
% 	W   Input argument
%
% OUTPUT
% 	N   1/0 if W is/isn't a mapping object
%
% DESCRIPTION
% True (1) if W is a mapping object and false (0), otherwise.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% ISDATASET, ISFEATIM, ISDATAIM

% $Id: ismapping.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function n = ismapping(w)
		n = isa(w,'prmapping');

	if (nargout == 0) & (n == 0)
		error([newline 'Mapping expected.'])
	end
return;
