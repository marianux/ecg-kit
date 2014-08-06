%PRDATAFILES Checks availability of a PRTools datafile (PRTools5 version!)
%
%   PRDATAFILES
%
% Checks the availability of the PRDATAFILES directory, downloads the
% Contents file and m-files if necessary and adds it to the search path. 
% Lists Contents file.
%
%   PRDATAFILES RENEW
%
% Replace PRDATAFILES m-files by its most recent version.
%
%   PRDATAFILES ALL
%
% Download and save all data related to the m-files (very time consuming!)
%
%		PRDATAFILES(DFILE)
%
% Checks the availability of the particular datafile DFILE. DFILE should be
% the name of the m-file. If it does not exist in the 'prdatafiles'
% directory an attempt is made to download it from the PRTools web site.
%
%   PRDATAFILES(DFILEDIR,SIZE,URL)
%
% This command should be used inside a PRDATAFILES m-file. It checks the 
% availability of the particular datafile directory DFILEDIR and downloads 
% it if needed. SIZE is the size of the datafile in Mbyte, just used to 
% inform the user. In URL the web location may be supplied. Default is 
% http://prtools.org/prdatafiles/DFILEDIR.zip
%
% All downloading is done interactively and should be approved by the user.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATAFILES, PRDATASETS, PRDOWNLOAD

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function prdatafiles(dfile,siz,url)

persistent ask;
if isempty(ask), ask = true; end 

if nargin < 1, dfile = []; end
if nargin < 2, siz = []; end
if nargin < 3 | isempty(url)
	url = ['http://prtools.org/prdatafiles/' dfile '.zip']; 
end

if exist('highway','file') ~= 2
  prtoolsdir = fileparts(which(mfilename));
  toolsdir = fileparts(prtoolsdir);
  if exist(fullfile(toolsdir,'prdatafiles/Contents.m'),'file') ~= 2
    path = input(['The directory prdatafiles is not found in the search path.' ... 
      newline 'If it exists, give the path, otherwise hit the return for an automatic download.' ...
      newline 'Path to prdatafiles: '],'s');
    if ~isempty(path)
      addpath(path);
      feval(mfilename,dfile,siz);
      return
    else
      % Load all m-files from prdatafiles5 !!!
      [ss,dirname] = prdownload('http://prtools.org/prdatafiles5/prdatafiles.zip','prdatafiles');
      addpath(dirname)
    end
  else
    addpath(fullfile(toolsdir,'prdatafiles'));
  end
end

if isempty(dfile) % just list Contents file
	
	help('prdatafiles/Contents')
	
elseif ~isempty(dfile) & nargin == 1 % check / load m-file
	% this just loads the m-file in case it does not exist and updates the
	% Contents file
	if strcmpi(dfile,'renew')
		if exist('prdatafiles/Contents.m','file') ~= 2
			% no prdatafiles in the path, just start
			feval(mfilename);
		else
			dirname = fileparts(which('prdatafiles/Contents'));
			prdownload('http://prtools.org/prdatafiles5/prdatafiles.zip',dirname);
    end
  elseif strcmpi(dfile,'all')
		if exist('prdatafiles/Contents','file') ~= 2
			% no prdatafiles in the path, just start
			feval(mfilename);
    end
    % load all data without asking
    ask = false;
    tooldir = fileparts(which('prdatafiles/Contents'));
    files = dir([tooldir '/*.m']);
    files = char({files(:).name});
    L = strmatch('Contents',files); % no data for Contents
    L = [L; strmatch('pr',files)];  % no data for support routines
    files(L,:) = [];
    for j=1:size(files,1)
      cmd = deblank(files(j,:));
      disp([newline cmd])
      feval(cmd(1:end-2));
    end
    ask = true;
    
	elseif exist(['prdatafiles/' dfile],'file') ~= 2 
    % load m-file
		prdownload(['http://prtools.org/prdatafiles5/' dfile '.m'],dirname);
		prdownload('http://prtools.org/prdatafiles5/Contents.m',dirname);
    feval(dfile);   % takes care that data is available as well
	end
	
else   % dfile is now the name of the datafile directory
	     %It might be different from the m-file, so we cann't check it.
	rootdir = fileparts(which('prdatafiles/Contents'));
	if exist(fullfile(rootdir,dfile),'dir') ~= 7
    if ask
      csiz = ['(' num2str(siz) ' MB)'];
      q = input(['Datafile is not available, OK to download ' csiz ' [y]/n ?'],'s');
      if ~isempty(q) & ~strcmp(q,'y')
        error('Datafile not found')
      end
    end
		prdownload(url,rootdir,siz);
		disp(['Datafile ' dfile ' ready for use'])
	end
end
