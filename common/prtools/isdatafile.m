%ISDATAFILE Test whether the argument is a datafile
%
% 	N = ISDATAFILE(A);
%
% INPUT
%		A	 Input argument
%
% OUTPUT
%		N  1/0 if A is/isn't a datafile
%
% DESCRIPTION
% The function ISDATAFILE test if A is a datafile object.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% ISMAPPING, ISDATAIM, ISFEATIM 

function n = isdatafile(a)
			
	n = isa(a,'prdatafile');
	if (nargout == 0) & (n == 0)
		error([newline 'Datafile expected.'])
	end
return;
