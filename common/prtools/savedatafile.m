%SAVEDATAFILE Save datafile
%
%		B = SAVEDATAFILE(A,FEATSIZE,NAME,NBITS,FILESIZE)
%
% INPUT
%   A          Datafile, or cell array with datafiles and/or datasets
%   FEATSIZE   Feature size, i.e. image size of a single object in B
%              Default: as it is.
%   NAME       Desired name of directory
%   NBITS      # of bits in case of rescaling (8,16 or 32)
%              Default: no rescaling
%   FILESIZE   # of elements stored in a single file
%              Default 10000000.
%
% OUTPUT
%    B         New datafile
%
% DESCRIPTION
% The datafile A is completed by desired preprocessing and postprocessing.
% It should therefore be convertable to a dataset. It is saved in the 
% directory NAME and B refers to it. If desired the data and target fields
% are compressed (after appropriate rescaling) to NBITS unsigned integers. 
% The stored datafile can be retrieved by
% 
%    B = PRDATAFILE(NAME)
%
% B is a 'mature' datafile, i.e. a dataset distributed over a number of
% files with maximum size FILESIZE. This has only advantages over a 'raw' 
% datafile defined for a directory of images in case of substantial pre- 
% and postprocessing, due to the overhead of the dataset construct of B. 
% FEATSIZE can be used to reshape the size of object (e.g. from 256 to
% [16 16])
%
% If A is cell array with datafiles and/or datasets, they are first
% horizontally concatenated before the datafile is written. The first
% element in A should be a datafile.
%
% The difference with CREATEDATAFILE is that SAVEDATAFILE assumes that
% the datafile after completion by preprocessing can be converted into a
% dataset and the data is stored as such and a sompact as possible.
% CREATEDATAFILE saves every object in a separate file and is thereby
% useful in case preprocessing does not yield a proper dataset with the
% same number of features for every object.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATAFILES, DATASETS, CREATEDATAFILE

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function c = savedatafile(dfile_in,featsize,dirr,nbits,filesize)
    
  if iscell(dfile_in)
    b = dfile_in{1};
    k = getsize(prdataset(b(1,:)),2);
    m = getsize(b,1);
    for j=2:length(dfile_in(:))
      if size(dfile_in{j},1) ~= m
        error('Datafiles or datasets to be combined should have equal number of objects')
      end
      k = k+size(dfile_in{j},2);
    end
  else
    b = dfile_in;
    k = getsize(prdataset(b(1,:)),2);
    m = getsize(b,1);
  end
     
	isvaldfile(b);
  c = b;

	if nargin < 5                   % maximum nr of elements in datafile
		filesize = 1000000;
	end
	if nargin < 4 | isempty(nbits)  % default: no conversion
		nbits = 0;
	end
  
	if nargin < 3 | isempty(dirr)   % invent name if not supplied
		name = getname(b);
		if ~isempty(name) & length(name) > 10
			name = [deblank(name(1:10)) '_'];
		else
			name = 'PRDatafile_';
		end
		dirr = [name num2str(round(rand*100000))]; % make it unique
		dirr = fullfile(pwd,dirr);  %DXD store it in the working dir
	end
	if all(nbits ~= [0 1 8 16 32])
		error('Illegal number of bits supplied for conversion')
	end
	
	if nargin < 2, featsize = []; end
	                      % reorder datafile for faster retrieval
  if isdatafile(b) & ~iscell(dfile_in)
    file_index = getident(b.prdataset,'file_index');
	  [jj,R] = sort(file_index(:,1)+file_index(:,2)/(max(file_index(:,2),[],1)+1));
    b = b(R,:);
  end
	%[m,k] = size(b);
	nobj= floor(filesize/k);   % # of objects per file
	if nobj<1
		error('Filesize is too small. A single object does not even fit in!');
	end
	nfiles = ceil(m/nobj);     % # of files needed
	if exist(dirr) == 7        % if datafile already exists
		act_dir = pwd;
		cd(dirr);                % go to it
		fclose('all');
		delete *.mat             % make it empty
		delete file*
		cd(act_dir)              % return to original directory
	else                       % otherwise create it
		[pp,ff] = fileparts(dirr);
    if isempty(pp)
      pp = pwd;
    end
    mkdir(pp,ff)
	end
  
	obj1 = 1;                  % initialise start object
	obj2 = min(obj1+nobj-1,m); % initialise end object
	file_index = zeros(m,2);
	s = sprintf('Saving %4i datafiles: ',nfiles);
	prwaitbar(nfiles,s);
	spaces = '        ';   
                             % initialise structure with file information
	cfiles = struct('name',cell(1,nfiles),'nbits',[], ...
        'offsetd',[],'scaled',[],'offsett',[],'scalet',[]);
	for j=1:nfiles             % run over all files to be created
		prwaitbar(nfiles,j,[s int2str(j)]);
		a = prdataset(b(obj1:obj2,:)); % convert objects for a single file into dataset
		if iscell(dfile_in)
			for i=2:length(dfile_in(:))
				a = [a prdataset(dfile_in{i}(obj1:obj2,:))];
			end
		end
		if ~isempty(featsize)
			if prod(featsize) ~= prod(getfeatsize(a))
				error('Desired feature size does not match size of data')
			else
				a = setfeatsize(a,featsize);
			end
		end
		lendata = size(a,2);
		[a,scaled,offsetd,scalet,offsett,prec] = convert(a,nbits); % convert precision
		filej = ['file',num2str(j)];
		fname = fullfile(dirr,filej);
   		fid = fopen(fname,'w');
		if fid<0
			error('Unable to write file %s.',fname);
		end
		fwrite(fid,[a.data a.targets]',prec);  % store
		fclose(fid);
    
		cfiles(j).name    = filej;      % name of the file
		cfiles(j).prec    = prec;       % precision
		cfiles(j).scaled  = scaled;     % scaling for datafield
		cfiles(j).offsetd = offsetd;    % offset for datafield
		cfiles(j).scalet  = scalet;     % scaling for target field
		cfiles(j).offsett = offsett;    % offset for target field
		cfiles(j).sized   = a.featsize; % feature (image) size
		cfiles(j).sizet   = size(a.targets,2); % size target field
                                    % create file_index
		file_index(obj1:obj2,1) = repmat(j,obj2-obj1+1,1);
		file_index(obj1:obj2,2) = [1:obj2-obj1+1]';
		obj1 = obj2+1;                  % update start object
		obj2 = min(obj1+nobj-1,m);      % update end object
	end
	prwaitbar(0);
                                    % finalise datafile definition
	c.files = cfiles;           % file information
	c.type = 'mature';                % mature datafile
	preproc.preproc = [];
	preproc.pars = {};
	c.preproc = preproc;              % no preprocessing
	c.postproc = prmapping([]);         % default postprocessing
	featsize = getfeatsize(c);        % copy feature (image) size
	c.postproc = setsize_in(c.postproc,featsize);  % set input and ...
	c.postproc = setsize_out(c.postproc,featsize); % output size postprocessing
	c.prdataset = setfeatsize(c.prdataset,getfeatsize(c)); % set size datafile
	c.prdataset = setident(c.prdataset,file_index,'file_index'); % store file_index
	%DXD give it the name of the directory, not the complete path:
	%c.prdataset = setname(c.prdataset,dirr);
	[pname,dname] = fileparts(dirr);
	c.prdataset = setname(c.prdataset,dname);
  if isempty(pname), pname = pwd; end
  c.rootpath = fullfile(pname,dname); 
  
	%DXD this path/filename is *very* strange!
	%save(fullfile(pwd,dirr,dirr),'c'); % store the datafile as mat-file for fast retrieval
	%DXD I changed it into this:
	%save(fullfile(dirr,filej),'c'); % store the datafile as mat-file for
	%fast retrieval
	save(fullfile(dirr,dname),'c'); % store the datafile as mat-file for fast retrieval
  
	%closemess([],len1+len2);
	
return

function [a,scaled,offsetd,scalet,offsett,prec] = convert(a,nbits)
                    % handle precision data and target fields
  if nbits == 0     % no conversion
    scaled = [];
    offsetd = [];
    scalet = [];
    offsett = [];
    prec = 'double';
    return; 
  end
                    % compute data field conversion pars
  data = +a;
  maxd = max(max(data));
	mind = min(min(data));
	scaled = (2^nbits)/(maxd-mind);
	offsetd = -mind;
	data = (data+offsetd)*scaled;
	                  % compute target field conversion pars
	targ = gettargets(a);
	if isempty(targ)  % no targets set
		scalet = [];
		offsett = [];
	else
		maxt = max(max(targ));
		mint = min(min(targ));
		scalet = (2^nbits)/(maxt-mint);
		offsett = -mint;
		targ = (targ+offsett)*scalet;
	end
	                   % execute conversion
	switch nbits
		case 8
			data = uint8(data);
			targ = uint8(targ);
      prec = 'uint8';
		case 16
			data = uint16(data);
			targ = uint16(targ);
      prec = 'uint16';
		case 32
			data = uint32(data);
			targ = uint32(targ);
      prec = 'uint32';
		otherwise
			error('Desired conversion type not found')
	end
                      % store results
  a.data = data;
  a.targets = targ;
return
	

