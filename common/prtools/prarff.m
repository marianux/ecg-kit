%PRARFF COnvert ARFF file into PRTools dataset
%
%		A = PRARFF(FILE)
%
% INPUT
%   FILE    ARFF file
%
% OUTPUT
%   A       Dataset in PRTools format
%
% DESCRIPTION
% ARFF files as used in WEKA are converted into PRTools format. In case
% they don't fit (non-numeric features, varying feature length) an error is
% generated.
%
% SEE ALSO
% DATASETS

function a = prarff(file)

if nargin < 1 | exist(file) ~= 2
	error('file not found');
end

t = txtread(file);
c = cell2str(t);

k = 0;
nodata = 1;
for j=1:length(c);
	if nodata
		[s,u] = strtok(c{j});
		if strcmp(s,'@relation')
			name = strtrim(u);
		end
		if strcmp(s,'@attribute')
			k = k+1;
			u = strrep(u,'''','');
			[featlab{k},u] = strtok(u);
 			if strcmp(featlab{k},'class')
				featlab(k) = [];
 				k = k-1;
%         skip as we determine lablist from data field
% 				u = strrep(u,'{','{''');
% 				u = strrep(u,'}','''}');
% 				u = strrep(u,',',''',''');
% 				eval(['lablist = ' u]);
			else
				if ~any(strcmp(strtrim(u),{'numeric','integer','real'}))
					warning(['Non-numeric attributes are not supported ' ...
						file ' ' u]);
					a = []; 
					return
				end
			end
		end
		if strcmp(s,'@data')
			nodata = 0;
			form = [repmat('%e,',1,k) '%s'];
			a = zeros(length(c)-j,k);
			m = 0;
		end
	else
		if length(find(c{j}==',')) == k
			m = m+1;
			x = sscanf(c{j},form);
			a(m,:) = x(1:k);
			labels{m} = char(x(k+1:end))';
		elseif length(find(c{j}==',')) == k-1 % unlabeled?
			m = m+1;
			x = sscanf(c{j},form);
			a(m,:) = x(1:k);
			labels{m} = '';
		else
			error('Data size doesn''t match number of attributes')
		end
	end
end
a = prdataset(a,labels);
a = setfeatlab(a,featlab);
a = setname(a,name);

return
		
		
%STR2CELL String to cell conversion
%
%		C = STR2CELL(S)
%
% INPUT
%   S    String
%
% OUTPUT
%   A    Cell array
%
% DESCRIPTION
% The string S is broken into a set of strings, one for each line. Each of
% them is place into a different element of the cell araay C

function c = cell2str(s)

if nargin < 1 | ~ischar(s)
	error('No input string found')
end

s = strrep([s char(10)],char([10 13]),char(10));
s = strrep(s,char([13 10]),char(10));
s = strrep(s,char([10 10]),char(10));
s = strrep(s,char(13),char(10));
n = [0 strfind(s,char(10))];

c = cell(length(n-1),1);
for j=1:length(n)-1
	c{j} = s(n(j)+1:n(j+1)-1);
end
if isempty(c{end})
	c(end) = [];
end

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
			
		
	

% 
