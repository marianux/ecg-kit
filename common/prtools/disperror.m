%DISPERROR Display error matrix with information on classifiers and datasets
%
%	DISPERROR(DATA,CLASSF,ERROR,STD,FID)
% 
% INPUT
%   DATA     Cell array of M datasets or dataset names (strings)
%   CLASSF   Cell array of N mappings or mapping names (strings)
%   ERROR    M*N matrix of (average) error estimates 
%   STD      M*N matrix of standard devations on ERROR (optional)
%   FID      File in which results are written (default: 1)
% OUTPUT
%
% DESCRIPTION
% Displays the matrix ERROR matrix with error estimates for N
% classifiers related to M datasets. This routine is called by TESTC
% and CROSVALL to display results.
%
% EXAMPLE
% testsets  = {gendath gendatb gendatd(100,5)}
% trainsets = {gendath gendatb gendatd(100,5)}
% classifiers = {nmc fisherc qdc svc}
% testc(testsets,prmap(trainsets,classifiers))
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, TESTC, PRCROSSVAL

% $Id: disperror.m,v 1.3 2007/06/05 12:43:35 duin Exp $

function disperror (data,classf,err,stdev,fid)
	
	if nargin < 5, fid = 1; end
	% Check arguments.
	if (nargin > 3) & (any(size(err) ~= size(stdev)))
		error('size of matrix with standard deviations should match matrix with errors')
	end
	if (~iscell(classf)) | (~isstr(classf{1}) & ~ismapping(classf{1}))
		error('cell array of mappings or mapping names expected')
	end
	if (~iscell(data)) | (~isstr(data{1}) & ~isdataset(data{1}))
		error('cell array of datasets or datasets names expected')
	end

	[m,n] = size(err);
	if (length(data) ~= m)
		error('size of dataset cell array should equal number of rows in error matrix');
	end

	if (length(classf) ~= n)
		error('size of classifier cell array should equal number of columns in error matrix');
	end

	% If datasets are supplied, extract their names.
	for j = 1:m
		if (isdataset(data{j}))
			data{j} = getname(data{j});
		end
	end

	% If classifiers are supplied, extract their names.
	for j = 1:n
		if (ismapping(classf{j}))
			classf{j} = getname(classf{j});
		end
  end

  if n >=  m
    
    if m == 1
      fprintf(fid, ' %s \n\n',data{1});
    else
      fprintf(fid,'\n');
      for j = 1:m
        fprintf(fid,'\n  data_%i : %20s',j,data{j});
      end
      fprintf(fid,'\n\n                      ');
      for j = 1:m
        fprintf(fid,'  data_%i',j);
      end
      fprintf(fid,'\n\n');
    end

    for i = 1:n
      fprintf(fid,'  %-22s',classf{i});
      fprintf(fid,'  %5.3f',err(:,i)');
      if (nargin > 3)
        fprintf(fid,' (%5.3f)',stdev(:,i)');
        fprintf(fid,'\n');
      end
      fprintf(fid,'\n');
    end
    
  else
      
    if (n == 1)
      fprintf(fid,' %s \n\n',classf{1});
    else
      fprintf(fid,'\n');
      for i = 1:n
        fprintf(fid,'\n  clsf_%i : %s',i,classf{i});
      end
      fprintf(fid,'\n\n                      ');
      for i = 1:n
        fprintf(fid,'  clsf_%i',i);
      end
      fprintf(fid,'\n\n');
    end

    for j = 1:m
      fprintf(fid,'  %s',data{j});
      fprintf(fid,' %7.3f',err(j,:));
      if (nargin > 3)
        fprintf(fid,'\n                      ');
        fprintf(fid,' %7.3f',stdev(j,:));
        fprintf(fid,'\n');
      end
      fprintf(fid,'\n');
    end
    
  end
  fprintf(fid,'\n');
	
return
