%NODATAFILE Error return in case of datafile
%
%		NODATAFILE
%
% Error message
%
%		B = NODATAFILE(A)
%		B = A*NODATAFILE
%
% Error message in case A is a datafile, otherwise B = A

function  a = nodatafile(a)

	if (nargin == 0 & nargout == 0) | (nargin == 1 & isdatafile(a) & nargout == 0)
		error('prtools:nodatafile','Command not implemented for datafiles');
	elseif nargin == 0
		a = prmapping(mfilename,'fixed');
	else
		;
	end
	
return;
 