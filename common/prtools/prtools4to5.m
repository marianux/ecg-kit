%PRTOOLS4TO5
%
%   PRTOOLS4TO5(DIRIN,DIROUT)
%
% Copies all files in the directory DIRIN into DIROUT.
% Convert all occurances in DIROUT of calls to
%
%    dataset, datafile, mapping, map, pca and crossval  
%
% to
%
%    prdataset, prdatafile, prmapping, prmap, pcam and prcrossval
%
% This routine may not be full proof. Result should be checked.

% Copyright: R.P.W. Duin, r.p.w.duin37steps.com

function prtools4to5(dirin,dirout)

if exist(dirout,'dir') == 7
  [suc,mess]=rmdir(dirout,'s');
  disp(mess)
end
mkdir(dirout);

copydir(dirin,dirout);
replacem(dirout,' dataset(',' prdataset(');
replacem(dirout,' dataset;',' prdataset;');
replacem(dirout,'+dataset(','+prdataset(');
replacem(dirout,'=dataset(','=prdataset(');
replacem(dirout,',dataset)',',prdataset)');
replacem(dirout,',dataset(',',prdataset(');
replacem(dirout,'(dataset(','(prdataset(');
replacem(dirout,'.dataset','.prdataset');
replacem(dirout,'''dataset''','''prdataset''');
replacem(dirout,' DATASET(',' PRDATASET(');
replacem(dirout,' DATASET,',' PRDATASET,');
replacem(dirout,',DATASET,',',PRDATASET,');
replacem(dirout,' DATASET.',' PRDATASET.');
replacem(dirout,'%DATASET ','%PRDATASET ');

replacem(dirout,' mapping(',' prmapping(');
replacem(dirout,' mapping;',' prmapping;');
replacem(dirout,'+mapping(','+prmapping(');
replacem(dirout,'*mapping(','*prmapping(');
replacem(dirout,'=mapping(','=prmapping(');
replacem(dirout,',mapping)',',prmapping)');
replacem(dirout,',mapping(',',prmapping(');
replacem(dirout,'.mapping ','.prmapping ');
replacem(dirout,'.mapping=','.prmapping=');
replacem(dirout,'.mapping,','.prmapping,');
replacem(dirout,'''mapping''','''prmapping''');
replacem(dirout,' MAPPING(',' PRMAPPING(');
replacem(dirout,' MAPPING,',' PRMAPPING,');
replacem(dirout,',MAPPING,',',PRMAPPING,');
replacem(dirout,' MAPPING.',' PRMAPPING.');
replacem(dirout,'%MAPPING ','%PRMAPPING ');

replacem(dirout,' datafile(',' prdatafile(');
replacem(dirout,' datafile;',' prdatafile;');
replacem(dirout,'+datafile(','+prdatafile(');
replacem(dirout,'=datafile(','=prdatafile(');
replacem(dirout,',datafile)',',prdatafile)');
replacem(dirout,'.datafile ','.prdatafile ');
replacem(dirout,'.datafile=','.prdatafile=');
replacem(dirout,'.datafile,','.prdatafile,');
replacem(dirout,'''datafile''','''prdatafile''');
replacem(dirout,' DATAFILE(',' PRDATAFILE(');
replacem(dirout,' DATAFILE,',' PRDATAFILE,');
replacem(dirout,',DATAFILE,',',PRDATAFILE,');
replacem(dirout,' DATAFILE.',' PRDATAFILE.');
replacem(dirout,'%DATAFILE ','%PRDATAFILE ');

replacem(dirout,' map(',' prmap(');
replacem(dirout,'+map(','+prmap(');
replacem(dirout,'*map(','*prmap(');
replacem(dirout,'=map(','=prmap(');
replacem(dirout,',map(',',prmap(');
replacem(dirout,' MAP(',' PRMAP(');
replacem(dirout,' MAP,',' PRMAP,');
replacem(dirout,',MAP,',',PRMAP,');
replacem(dirout,' MAP.',' PRMAP.');
replacem(dirout,'%MAP ','%PRMAP ');

replacem(dirout,' pca(',' pcam(');
replacem(dirout,'''pca''','''pcam''');
replacem(dirout,'+pca(','+pcam(');
replacem(dirout,'*pca(','*pcam(');
replacem(dirout,'=pca(','=pcam(');
replacem(dirout,',pca(',',pcam(');
replacem(dirout,';pca(',';pcam(');
replacem(dirout,' PCA(',' PCAM(');
replacem(dirout,' PCA,',' PCAM,');
replacem(dirout,',PCA,',',PCAM,');
replacem(dirout,',PCA ',',PCAM ');
replacem(dirout,' PCA.',' PCAM.');
replacem(dirout,'%PCA ','%PCAM ');

replacem(dirout,' crossval(',' prcrossval(');
replacem(dirout,'=crossval(','=prcrossval(');
replacem(dirout,',crossval(',',prcrossval(');
replacem(dirout,';crossval(',';prcrossval(');
replacem(dirout,' CROSSVAL(',' PRCROSSVAL(');
replacem(dirout,' CROSSVAL,',' PRCROSSVAL,');
replacem(dirout,',CROSSVAL,',',PRCROSSVAL,');
replacem(dirout,',CROSSVAL ',',PRCROSSVAL ');
replacem(dirout,' CROSSVAL.',' PRCROSSVAL.');
replacem(dirout,'%CROSSVAL ','%PRCROSSVAL ');
return

function replacem(dir,s1,s2)
% replace in dir and all subdirs string s1 by s2. m-files only
[subdirs,files] = dirnames(dir);
if ~isempty(subdirs)
  for j=1:size(subdirs,1)
    subname = deblank(subdirs(j,:));
    replacem(fullfile(dir,subname),s1,s2);
  end
end
if ~isempty(files)
  for j=1:size(files,1)
    filename = deblank(files(j,:));
    [~,~,ext] = fileparts(filename);
    if strcmp(ext,'.m')
      replacef(fullfile(dir,filename),s1,s2);
    end
  end
end

function replacef(file,s1,s2)
%replace in text file all occurances of string s1 by s2
r = readf(file);
n = grep(r,s1);
if ~isempty(n)
  disp(['repl ' file])
  c = listn(r);               % lines in cell array
  c(n) = strrep(c(n),s1,s2);  % convert the lines of interest
  r = [c{:}];                 % back to a single string
  writf(file,r);
end
return

%COPYDIR Copy all files from dir to dir
function copydir(dir1,dir2)
if exist(dir2,'dir') ~= 7
	mkdir(dir2);
end

[subdirs,files] = dirnames(dir1);
if ~isempty(subdirs)
  for j=1:size(subdirs,1)
    subname = deblank(subdirs(j,:));
    copydir(fullfile(dir1,subname),fullfile(dir2,subname));
  end
end
if ~isempty(files)
  for j=1:size(files,1)
    filename = deblank(files(j,:));
    disp(['copy ' filename])
    copyfile(fullfile(dir1,filename),fullfile(dir2,filename),'f');
  end
end
    
%DIRNAMES
%
%	[SUBDIRS,FILES]= DIRNAMES(DIR)
%
%Get names of all subdirs and files in direcotory DIR.
%Hidden files, . and .. neglected

function [subdirs,files] = dirnames(dirname)

allnames = dir(dirname);
n = length(allnames);
subdirs = cell(1,n);
files = cell(1,n);
ns = 0;
nf = 0;

for j=1:n
	if allnames(j).name(1) ~= '.' % skip hidden files
		if allnames(j).isdir
			ns = ns+1;
			subdirs{ns} = allnames(j).name;
		else
			nf = nf+1;
			files{nf} = allnames(j).name;
		end
	end
end

if ns > 0
	subdirs = char(subdirs(1:ns));
else 
	subdirs = [];
end
if nf > 0
	files = char(files(1:nf));
else
	files = [];
end

%GREP Get line specific lines
%
% [k,n] = grep(r,s)
% Get the numbers of all lines in the set of lines r 
% that contain s.
% n is the total number of lines.

function [k,z] = grep(r,s)
n = [0,find(r==newline)];
m = strfind(r,s);
[~,j] = sort([n,m]);
q = [0,j(1:length(j)-1)]-j;
k = j(q>0)-1;
z = length(n)-1; % # of lines
return

%NEWLINE The platform dependent newline character
%
% c = newline

function c = newline
	
	if strcmp(computer,'MAC2')
		c = char(13);
	elseif strcmp(computer,'PCWIN')
		c = char(10);
	else
		c = char(10);
	end
	
return

%WRITF Write file
%
% writf(file,r)
% Write file from string r

function writf(file,r)
fid = fopen(file,'w');
if fid < 0
   error('Cannot open file')
end
fprintf(fid,'%c',r);
fclose(fid);
return

%READF Readfile
%
% [r,n] = readf(file)
% Reads file into string r. The number of lines
% is returned in n.

function [r,n] = readf(file)
fid = fopen(deblank(file),'r');
if fid < 0
   error(['Cann''t open ' file])
end
r = fscanf(fid,'%c');
fclose(fid);
n = length(find(r==13));
if r(length(r)) ~= 13, n = n + 1; end
return

%LISTN List lines specified by their line number
%
% t = listn(r,n)
% Get the lines in r given by the line numbers in n.
% Default n, get all.
% t is a cell array
function t = listn(r,n)
k = [0,find(r==newline)];
if k(end)~=length(r)
  k = [k length(r)];
end
if nargin < 2, n = [1:length(k)-1]; end
lenr = length(r);
lenk = length(k);
lenn = length(n);
t = cell(1,length(n));
for j = 1:lenn
	if n(j) < lenk 
		m = k(n(j))+1:k(n(j)+1);
		if all(m <= lenr)
  		t{j} = r(m);
		end
	end
end
return