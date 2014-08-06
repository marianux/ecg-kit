%PRDOWNLOAD Download well defined data and toolboxes
%
%   STATUS = PRDOWNLOAD(URL,DIRNAME)
%
% INPUT
%   URL          String containing URL of file to be downloaded
%   DIRNAME      Final directory for download, created if necessary.
%                Default: directory of calling function
%
% DESCRIPTION
% The URL will be downloaded in directory DIRNAME (to be created if
% needed). The resulting file will be uncompressed in case of a zip-, gz- 
% or tar-file. 
%
% The main purpose of this routine is to download missing datafiles or
% datasets from PRDATAFILES and PRDATASETS.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% PRDATASETS, PRDATAFILES

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function [status,dirname] = prdownload(url,dirname,siz)

% create target directory if needed
if nargin < 3, siz = []; end
if nargin > 1 
  if exist(dirname,'dir') ~= 7
    success = mkdir(dirname);
    if ~success
    	maindir = input(['It is not possible to create a directory here.' newline ...
        'Please give another location: '],'s');
      dirname = fullfile(maindir,dirname);
      mkdir(dirname);
    end
  end
else
  dirname = callerdir;
end

% download
[dummy,ff,xx]= fileparts(url);
filename = [ff xx];
if isempty(siz)
  disp(['Downloading ' filename ' ....'])
else
  disp(['Downloading ' filename ' (' num2str(siz) ' MB) ....'])
end
rfilename = fullfile(dirname,filename);
if ~usejava('jvm') & isunix 
  disp('isunix')
  disp(['wget -q -O ' rfilename [' '] url])
	stat = unix(['wget -q -O ' rfilename [' '] url]);
	status = (stat == 0);
else
  if verLessThan('matlab','8')
    [f,status] = urlwrite(url,rfilename);
  else
    [f,status] = urlwrite(url,rfilename,'TimeOut',20);
  end
end
if status == 0
	error('Server unreachable, timed out or file not found')
end

decompress(rfilename);

function dirname = callerdir

ss = dbstack;
if length(ss) < 3
  % no caller, commandline call
  dirname = pwd;
else
  dirname = fileparts(ss(3).name);
end

function decompress(file)

[dirname,filename,ext] = fileparts(file);
if any(strcmp(ext,{'.zip','.gz','.tgz','.tar'}))
  disp('Decompression ....')
  if strcmp(ext,'.zip')
    if ~usejava('jvm') & isunix
      [stat,s] = unix(['unzip ' file ' -d ' dirname]);
    else
      unzip(file,dirname);
    end
  elseif strcmp(ext,'.gz')
    gunzip(file,dirname);
  elseif strcmp(ext,'.tar')| strcmp(ext,'.tgz') 
    untar(file,dirname);
  end
  delete(file);
  decompress(fullfile(dirname,filename));
end

