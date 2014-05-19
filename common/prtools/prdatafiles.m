%PRDATAFILES Checks availability of a PRTools datafile
%
%   PRDATAFILES
%
% Checks the availability of the 'prdatafiles' directory, downloads the
% Contents file and m-files if necessary and adds it to the search path. 
% Lists Contents file.
%
%   PRDATAFILES RENEW
%
% Replace PRDATASETS by its most recent version.
%
%   PRDATAFILES(DFILEDIR,SIZE,URL)
%
% This command should be used inside a PRDATAFILES m-file. It checks the 
% availability of the particular datafile directory DFILEDIR and downloads 
% it if needed. SIZE is the size of the datafile in Mbyte, just used to 
% inform the user. In URL the web location may be supplied. Default is 
% http://37steps.com/prdatafiles/DFILEDIR.zip
%
% All downloading is done interactively and should be approved by the user.
%
% SEE ALSO
% DATAFILES, PRDATASETS, PRDOWNLOAD

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function prdatafiles(dfile,siz,url)

if nargin < 1, dfile = []; end
if nargin < 2, siz = []; end
if nargin < 3 | isempty(url)
	url = ['http://37steps.com/prdatafiles/' dfile '.zip']; 
end

if exist('prdatafiles/Contents.m','file') ~= 2
	path = input(['The directory prdatafiles is not found in the search path.' ... 
		newline 'If it exists, give the path, otherwise hit the return for an automatic download.' ...
		newline 'Path to prdatafiles: '],'s');
	if ~isempty(path)
		addpath(path);
		feval(mfilename,dfile,siz);
		return
	else
		[ss,dirname] = prdownload('http://37steps.com/data/pfiles/prdatafiles/prdatafiles.zip','prdatafiles');
		addpath(dirname)
	end
end

if isempty(dfile) % just list Contents file
	
	help('prdatafiles/Contents')
	
elseif ~isempty(dfile) & nargin == 1 % check / load m-file
	% renew all m-files
	if strcmp(lower(dfile),'renew')
		if exist('prdatafiles/Contents.m','file') ~= 2
			% no prdatafiles in the path, just start
			feval(mfilename);
		else
			dirname = fileparts(which('prdatafiles/Contents'));
			prdownload('http://37steps.com/prdatafiles/prdatafiles.zip',dirname);
		end
	% this just loads the m-file in case it does not exist and updates the
	% Contents file
% 	elseif exist(['prdatafiles/' dfile],'file') ~= 2
% 		prdownload(['http://37steps.com/prdatafiles/' dfile '.m'],'prdatafiles');
% 		prdownload('http://37steps.com/prdatafiles/Contents.m','prdatafiles');
	end
	
else   % dfile is now the name of the datafile directory
	     %It might be different from the m-file, so we cann't check it.
	rootdir = fileparts(which('prdatafiles/Contents'));
	if exist(fullfile(rootdir,dfile),'dir') ~= 7
		siz = ['(' num2str(siz) ' MB)'];
		q = input(['Datafile is not available, OK to download ' siz ' [y]/n ?'],'s');
		if ~isempty(q) & ~strcmp(q,'y')
			error('Datafile not found')
		end
		prdownload(url,rootdir);
		disp(['Datafile ' dfile ' ready for use'])
	end
end
