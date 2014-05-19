%PRDOWNLOAD Download well defined data and toolboxes
%
%   STATUS = PRDOWNLOAD(URL,DIRNAME)
%
% INPUT
%   URL          String containing URL of file to be downloaded
%   DIRNAME      Final directory for download, created if necessary
%
% DESCRIPTION
% The URL will be downloaded in directory DIRNAME (to be created if
% needed). The resulting file will be uncompressed in case of a zip-, gz- 
% or tar-file. 
%
% The main purpose of this routine is to download missing datafiles or
% datasets from PRDATAFILES and PRDATASETS.
%
% SEE ALSO
% PRDATASETS, PRDATAFILES

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [status,dirname] = prdownload(url,dirname)

% create target directory if needed
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
  dirname = callerdir
end

% download
[dummy,ff,xx]= fileparts(url);
filename = [ff xx];
disp(['Downloading ' filename ' ....'])
rfilename = fullfile(dirname,filename);
if ~usejava('jvm') & isunix 
	stat = unix(['wget -q -O ' rfilename [' '] url]);
	status = (stat == 0);
else
	[f,status] = urlwrite(url,rfilename);
end
if status == 0
	error('Server unreachable or file not found')
end

% assume file is created, uncompress if needed
% delete compressed file
if strcmp(xx,'.zip')
	disp('Decompression ....')
	if ~usejava('jvm') & isunix
		[stat,s] = unix(['unzip ' rfilename ' -d ' dirname]);
	else
		unzip(rfilename,dirname);
	end
	delete(rfilename);
elseif strcmp(xx,'.gz')
	disp('Decompression ....')
	gunzip(filename,dirname);
	delete(filename);
elseif strcmp(xx,'.tar')| strcmp(xx,'.tgz') | strcmp(xx,'.tar.gz')
	disp('Decompression ....')
	untar(filename,dirname);
	delete(filename);
end
dirname = fullfile(pwd,dirname);
disp('Ready')

function dirname = callerdir

[ss ,i] = dbstack;
if length(ss) < 3
	fil = [];
else
	fil = ss(3).name;
end
dirname = fileparts(fil);



