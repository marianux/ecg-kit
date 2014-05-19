%TXTREAD Read text file
% 
% 	A = TXTREAD(FILENAME,N,START)
% 
% INPUT
%   FILENAME  Name of delimited ASCII file
%   N         Number of elements to be read (default all)
%   START     First element to be read (default 1)
%  
% OUTPUT
%   A         String
% 
% DESCRIPTION
% Reads the total file as text string into A

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function a = txtread(file,n,nstart)

	if nargin < 3, nstart = 1; end
	if nargin < 2 | isempty(n), n = inf; end

	fid = fopen(file);
	if (fid < 0)
		error('Error in opening file.')
	end
	a = fscanf(fid,'%c',n);
	fclose(fid);
	
return