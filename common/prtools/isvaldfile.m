%ISVALDFILE Test whether the argument is a valid datafile or dataset
%
% 	N = ISVALDFILE(A);
% 	N = ISVALDFILE(A,M);
% 	N = ISVALDFILE(A,M,C);
%
% INPUT
%   A     Input argument, to be tested on datafile or dataset
%   M     Minimum number of objects per class in A
%   C     Minimum number of classes in A
%
% OUTPUT
%   N  1/0 if A is / isn't a valid datafile or dataset
%
% DESCRIPTION
% The function ISVALDFILE tests if A is a datafile or dataset that has
% at least C classes. It is an extension of ISVALDSET and can be used it 
% when datasets as well as datafiles are allowed.
% For datafiles(sets) with soft or targets labels it is tested whether A 
% has at least M objects.
%
% SEE ALSO
% ISDATASET, ISMAPPING, ISDATAIM, ISFEATIM, ISCOMDSET, ISVALDSET

function n = isvaldfile(a,m,c)
		
	if nargin < 3 | isempty(c), c = 0; end
	if nargin < 2 | isempty(m), m = 0; end
	
	if nargout == 0
		if ~isdatafile(a)
			isdataset(a);
		end
	else
		n = 1;
		if ~isdatafile(a) & ~isdataset(a)
			n = 0;
			return
		end
	end
	
	if c == 1 & getsize(a,3) == 0
		;  % accept unlabeled dataset as having one class
	elseif c > 0 & getsize(a,3) == 0
		
		if nargout == 0
  		error([newline 'Labeled prdatafile(set) expected'])
		else
			n = 0;
		end
		
	elseif getsize(a,3) < c
		
		if nargout == 0
			error([newline 'Datafile(set) should have at least ' num2str(c) ' classes'])
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
		error([newline 'Datafile(set) should have at least ' num2str(m) ' objects'])
	end
	
return
