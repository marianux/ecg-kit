%ISCOMDSET Test whether datasets are compatible
%
% 	N = ISCOMDSET(A,B,CLAS);
%
% INPUT
%   A    Input argument, to be tested on dataset
%   B    Input argument, to be tested on compatibility with A
%   CLAS 1/0, test on equal classes (1) or don't test (0)
%        (optional; default 1)
%
% OUTPUT
%   N    1/0 if A and B are / are not compatible datasets
%
% DESCRIPTION
% The function ISCOMDSET tests whether A and B are compatible
% datasets, i.e. have the same features and the same classes.
%
% SEE ALSO
% ISDATASET, ISMAPPING, ISDATAIM, ISFEATIM, ISVALDFILE, ISVALDSET

function n = iscomdset(a,b,clas)
	
	if nargin < 3, clas = 1; end
	
	if isempty(b) % return of second dataset empty (i.e. not supplied)
		return
	end
	
	if nargout == 0
		isdataset(a);
		isdataset(b);
	else
		n = isdataset(a) & isdataset(b);
	end

	featlaba = setstr(getfeatlab(a));
	featlabb = setstr(getfeatlab(b));

	nf = strcmp(featlaba,featlabb);
		
	if ~nf
		if nargout == 0
			error([newline 'Datasets for training and testing/tuning should' newline ...
					'have the same features in the same order.'])
		else
			n = 0;
		end
	end
	
	if clas
		lablista = setstr(getlablist(a));
		lablistb = setstr(getlablist(b));

		no = strcmp(lablista,lablistb);
		
		if ~no
			if nargout == 0
				error([newline 'Datasets for training and testing/tuning should' newline ...
					'have the same class labels in the same order.'])
			else
				n = 0;
			end
		end
	end
	
return

