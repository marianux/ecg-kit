%DOUBLEM Datafile mapping for conversion to double
%
%	B = DOUBLEM(A)
%	B = A*DOUBLEM
%
% For datasets B = A, but the data inside A is converted to double.
% For datafiles B = A, but the command makes sure that the data in B is
% first converted to double when B is converted to a dataset.
% In all other cases A itself is converted to double.
% This is useful to convert INT8 and other formats inside a dataset.

function a = doublem(a)
		
	if nargin < 1 || isempty(a)
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
