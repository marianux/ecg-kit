%ISVALDSET Test whether the argument is a valid dataset
%
% 	N = ISVALDSET(A);
% 	N = ISVALDSET(A,M);
% 	N = ISVALDSET(A,M,C);
%
% INPUT
%   A     Input argument, to be tested on dataset
%   M     Minimum number of objects per class in A
%   C     Minimum number of classes in A
%
% OUTPUT
%   N  1/0 if A is / isn't a valid dataset
%
% DESCRIPTION
% Test whether A is a proper dataset with at least C classes.
% For crisp datasets this function tests whether A is a dataset that has
% at least M objects per class and C classes.
% For datasets with soft or targets labels it is tested whether A has at
% least M objects.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% ISDATASET, ISMAPPING, ISDATAIM, ISFEATIM, ISCOMDSET

function n = isvaldset(a,m,c)
		
	if nargin < 3 | isempty(c), c = 0; end
	if nargin < 2 | isempty(m), m = 1; end
	
	if nargout == 0
		isdataset(a);
	else
		n = 1;
		if ~isdataset(a)
			n = 0;
			return
		end
	end
	
	if c > 0 & getsize(a,3) == 0
		
		if nargout == 0
  		error([newline 'Labeled dataset expected'])
		else
			n = 0;
		end
		
	elseif getsize(a,3) < c
		
		if nargout == 0
			error([newline 'Dataset should have at least ' num2str(c) ' classes'])
		else
			n = 0;
		end
		
	end

	if islabtype(a,'crisp') & any(classsizes(a) < m)
		
		if nargout == 0
			if m == 1
				
				error([newline 'Classes should have at least one object.' newline ...
						'Remove empty classes by A = remclass(A).'])
			else
				cn = num2str(m);
				error([newline 'Classes should have at least ' cn ' objects.' ...
						newline 'Remove small classes by A = remclass(A,' cn ')'])
			end
		else
			n = 0;
		end
		
	end
	
	if islabtype(a,'soft','targets') & size(a,1) < m
		error([newline 'Dataset should have at least ' num2str(m) ' objects'])
	end

return
