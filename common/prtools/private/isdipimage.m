%ISDIPIMAGE Test whether dipimage is available
%
% 	N = ISDIPIMAGE;
%
% OUTPUT
%		N  1/0 if DIPIMAGE is available
%
% DESCRIPTION
% The function ISDIPIMAGE test whether DIPIMAGE is available.
% ISDIPIMAGE called without an output parameter generates an error
% if DIPIMAGE is not running.

function n = isdipimage
		
	n = exist('gaussf_adap_banana') == 2;
	
	if (nargout == 0) & (n == 0)
		error([newline 'DipImage is needed and not available'])
	end
	
return;
