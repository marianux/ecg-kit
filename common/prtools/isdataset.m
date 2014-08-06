%ISDATASET Test whether the argument is a dataset
%
% 	N = ISDATASET(A);
%
% INPUT
%		A	 Input argument
%
% OUTPUT
%		N  1/0 if A is/isn't a dataset
%
% DESCRIPTION
% The function ISDATASET test if A is a dataset object.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% ISMAPPING, ISDATAIM, ISFEATIM 

% $Id: isdataset.m,v 1.3 2007/03/22 08:54:59 duin Exp $

function n = isdataset(a)
			
	n = isa(a,'prdataset') & ~isa(a,'prdatafile');
	if (nargout == 0) & (n == 0)
		error([newline 'Dataset expected.'])
	end
return;
