%TESTDATASIZE of datafiles and convert to dataset
%
%	 B = TESTDATASIZE(A,STRING)
%  I = TESTDATASIZE(A,STRING,FALSE)
%  I = TESTDATASIZE(N)
%
% INPUT
%  A         datafile or dataset
%  STRING    'data' (default) or 'features' or 'objects'
%  N         Given data size to be tested
%
% OUTPUT
%  B         Converted dataset
%  I         TRUE:  conversion possible
%            FALSE: conversion not possible
%
% DESCRIPTION
% Depending on the value of PRMEMORY and the size of the datafile A, it is
% converted to a dataset, otherwise an error is generated.
% In case the third parameter is FALSE or the first is a scalar just a test
% is executed. In case of no output arguments an error is generated if
% conversion in impossible.
%
% The parameter STRING controls the type of comparison:
%
% 'data'          PROD(SIZE(A)) < PRMEMORY
% 'objects'       SIZE(A,1).^2  < PRMEMORY
% 'features'      SIZE(A,2).^2  < PRMEMORY
%
% SEE ALSO
% DATASETS, DATAFILES, PRMEMORY

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = testdatasize(a,type,flag)

	
	if nargin < 3,
		flag = 1;
	end
	
	if nargin < 2
		type = 'data';
  end
	
  b = true;
	if nargin == 1 & isdouble(a) & numel(a) == 1
    % estimated data size in a
    if a > prmemory
      if nargout == 0
        error(['Dataset too large for memory.' newline ...
          'Size (Mega elements) is ' int2str(round(a/1000000)) ', memory is ' int2str(prmemory/1000000) newline ...
          'Possible solutions:' newline ...
          '- decrease data size' newline ...
          '- increase PRMEMORY, see prmemory' newline ...
          '- consider batch processing, see setbatch']);
      else
        b = false;
      end
    end
    return
  end
  
	if isdataset(a) | isdouble(a)
		if flag
			b = a;
		end
		return
	end
	
  % Now we have a datafile
	a = setfeatsize(a,0); % featsize of datafiles is unreliable
	switch type
		case 'data'
			if prod(size(a)) > prmemory
				if flag
					error(['Dataset too large for memory.' newline ...
          'Size is ' int2str(prod(size(a))) ', memory is ' int2str(prmemory) newline ...
          'Possible solutions:' newline ...
          '- decrease data size' newline ...
          '- increase PRMEMORY, see prmemory' newline ...
          '- consider batch processing, see setbatch']);
				else
					b = false;
				end
			end
		case 'objects'
			if size(a,1).^2 > prmemory
				if flag
					error(['Number of objects too large for memory.' newline ...
          'Size is ' int2str(size(a,1).^2) ', memory is ' int2str(prmemory) newline ...
          'Possible solutions:' newline ...
          '- decrease data size' newline ...
          '- increase PRMEMORY, see prmemory' newline ...
          '- consider batch processing, see setbatch']);
				else
					b = false;
				end
			end
		case 'features'
			if size(a,2).^2 > prmemory
				if flag
					error(['Number of features too large for memory.' newline ...
          'Size is ' int2str(size(a,2).^2) ', memory is ' int2str(prmemory) newline ...
          'Possible solutions:' newline ...
          '- decrease data size' newline ...
          '- increase PRMEMORY, see prmemory' newline ...
          '- consider batch processing, see setbatch']);
				else
					b = false;
				end
			end
		otherwise
			error('Unknown test requested')
	end
	if nargout > 0
		if flag
			b = prdataset(a);
		end
	end
	
return
