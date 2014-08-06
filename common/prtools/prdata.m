%PRDATA Read data files
% 
% 	A = PRDATA(FILENAME,FLAG)
% 
% INPUT
%   FILENAME  Name of delimited ASCII file containing rows of data
%   FLAG      If not 0, first column is assumed to contain labels (default 1)
%  
% OUTPUT
%   A         Dataset
% 
% DESCRIPTION
% Reads data into the dataset A. The first word of each line is interpreted
% as label data. Each line is stored row-wise and interpreted as the feature 
% values of a single object. 
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, PRDATASET

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: prdata.m,v 1.3 2008/03/20 07:42:01 duin Exp $

function a = prdata(file,labels)

		% ASCII magic...
	if (strcmp(computer,'MAC2')), crlf = 13; else, crlf = 10; end

	if (nargin < 2)
		prwarning(4,'first column of file is assumed to contain labels');
	  labels = 1; 
	end

	% Open the file.

	fid = fopen(file);
	if (fid < 0)
		error('Error in opening file.')
	end

	% Read the data. First, find the number of items on the first line.

	s = fread(fid,inf,'uchar'); i = find(s==crlf);
	n = length(sscanf(setstr(s(1:i(1))),'%e'));

	fseek(fid,0,'bof');									% Return to the begin of the file.
	[a,num] = fscanf(fid,'%e',inf);			% Keep reading N objects per line.
	a = reshape(a,n,num/n)';						% Reshape to the correct data matrix.

	% Create the dataset, depending where the labels were stored.

	if (labels)
		lab=a(:,1); a(:,1)=[];
		a = prdataset(a,lab);
	else
		a = prdataset(a);
	end

  % And close the file.

	fclose(fid);

return
