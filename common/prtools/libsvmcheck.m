%LIBSVMCHECK  Check whether the LIBSVM toolbox is in the path

function libsvmcheck
	
	if exist('svmpredict','file') ~= 3 & exist('svmpredict.dll','file') ~= 2
		error([newline 'The LIBSVM package is not found. Download and install it from' ...
			newline 'http://www.csie.ntu.edu.tw/~cjlin/libsvm/'])
	else
		p = which('svmtrain');
		if ~isempty(strfind(p,'biolearn'))
			error('The bioinfo toolbox is in the way, please remove it from the path')
		end
	end
		