%PRDATASETS Checks availability of a PRTOOLS dataset
%
%   PRDATASETS
%
% Checks the availability of the PRDATASETS directory, downloads the
% Contents file and m-files if necessary and adds it to the search path. 
% Lists Contents file.
%
%   PRDATASET RENEW
%
% Replace PRDATASETS by its most recent version.
%
%   PRDATASETS ALL
%
% Download and save all datasets related to the m-files.
%
%		PRDATASETS(DSET)
%
% Checks the availability of the particular dataset DSET. DSET should be
% the name of the m-file. If it does not exist in the 'prdatasets'
% directory an attempt is made to download it from the PRTools web site.
%
%		PRDATASETS(DSET,SIZE,URL)
%
% This command should be used inside a PRDATASETS m-file. It checks the 
% availability of the particular dataset file and downloads it if needed. 
% SIZE is the size of the dataset in Mbyte, just used to inform the user.
% In URL the web location may be supplied.  Default is 
% http://37steps.com/prdatasets/DSET.mat
%
% All downloading is done interactively and should be approved by the user.
%
% SEE ALSO
% DATASETS, PRDATAFILES

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function prdatasets(dset,siz,url)

persistent ask;
if isempty(ask), ask = true; end 

if nargin < 3, url = []; end
if nargin > 0 & isempty(url)
	url = ['http://37steps.com/prdatasets/' dset '.mat']; 
end
if nargin < 2, siz = []; end
if nargin < 1, dset = []; end
dirname = fullfile(cd,'prdatasets');

if exist('prdatasets/Contents','file') ~= 6
	path = input(['The directory prdatasets is not found in the search path.' ... 
		newline 'If it exists, give the path, otherwise hit the return for an automatic download.' ...
		newline 'Path to prdatasets: '],'s');
	if ~isempty(path)
		addpath(path);
		feval(mfilename,dset,siz);
		return
	else
		%dirname = fullfile(cd,'prdatasets');
		[ss,dirname] = prdownload('http://37steps.com/data/pfiles/prdatasets/prdatasets.zip','prdatasets');
		addpath(dirname)
	end
end

if isempty(dset) % just list Contents file
	
	help('prdatasets/Contents')
	
elseif ~isempty(dset) & nargin == 1 % check / load m-file
	% this just loads the m-file in case it does not exist and updates the
	% Contents file
	if strcmp(lower(dset),'renew')
		if exist('prdatasets/Contents','file') ~= 6
			% no prdatasets in the path, just start
			feval(mfilename);
		else
			dirname = fileparts(which('prdatasets/Contents'));
			prdownload('http://37steps.com/prdatasets/prdatasets.zip',dirname);
    end
  elseif strcmp(lower(dset),'all')
		if exist('prdatasets/Contents','file') ~= 6
			% no prdatasets in the path, just start
			feval(mfilename);
    end
    % load all data without asking
    ask = false;
    tooldir = fileparts(which('prdatasets/Contents'));
    files = dir([tooldir '/*.m']);
    files = char({files(:).name});
    L = strmatch('Contents',files);
    L = [L; strmatch('pr_',files)];
    files(L,:) = [];
    for j=1:size(files,1)
      cmd = deblank(files(j,:));
      disp([newline cmd])
      feval(cmd(1:end-2));
    end
    ask = true;
    
    
	elseif exist(['prdatasets/' dset],'file') ~= 2
		prdownload(['http://37steps.com/prdatasets/' dset '.m'],dirname);
		prdownload('http://37steps.com/prdatasets/Contents.m',dirname);
	end
	
else   % now we load the m-file as well as the data given by the url
	
	% feval(mfilename,dset); % don't do this to allow for different mat-file
	% naming
	rootdir = fileparts(which('prdatasets/Contents'));
	[pp,ff,xx] = fileparts(url);
	if exist(fullfile(rootdir,[ff xx]),'file') ~= 2
    if ask
      siz = ['(' num2str(siz) ' MB)'];
      q = input(['Dataset is not available, OK to download ' siz ' [y]/n ?'],'s');
      if ~isempty(q) & ~strcmp(q,'y')
        error('Dataset not found')
      end
    end
		prdownload(url,rootdir);
		disp(['Dataset ' dset ' ready for use'])
	end
end
