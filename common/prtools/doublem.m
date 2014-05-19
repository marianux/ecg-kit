%DOUBLEM Datafile mapping for conversion to double
%
%	B = DOUBLEM(A)
%	B = A*DOUBLEM
%
% For datasets B = A, in all other cases A is converted to double.

function a = doublem(a)

		
	if nargin < 1 | isempty(a)
		a = prmapping(mfilename,'fixed');
		a = setname(a,'double');
	elseif isdataset(a)
		a = setdat(a,double(+a));
	elseif isdatafile(a)
		a = a*filtm([],'double');
	else
		a = double(a);
	end
	
return
