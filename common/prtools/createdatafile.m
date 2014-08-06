%CREATEDATAFILE Create datafile on disk
%
%		B = CREATEDATAFILE(A,DIR,ROOT,TYPE,CMD,FMT)
%
% INPUT
%   A          Datafile
%   DIR        Name of datafile, default: name of A
%   ROOT       Root directory in which datafile should be created
%              default: present directory
%   TYPE       Datafile type ('raw' (default) or 'cell')
%   CMD        Command for writing files, default: IMWRITE
%   FMT        Format, default 'bmp'
%
% OUTPUT
%    B         New datafile
%
% DESCRIPTION
% The existing datafile A is recreated on disk, one file per object,
% including all preprocessing. The IDENT.IDENT field and the class priors 
% of A are copied. All other addional information of A is lost.
% The order of objects in A and B may differ.
%
% Raw datafiles are stored as integer images, use TYPE = 'cell' to store
% floats.
%
% The difference between CREATEDATAFILE and SAVEDATAFILE is that the latter
% assumes that the datafile after completion by preprocessing can be converted 
% into a dataset and the data is stored as such and as compact as possible.
% CREATEDATAFILE saves every object in a separate file and is thereby
% useful in case preprocessing does not yield a proper dataset with the
% same number of features for every object.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATAFILES, DATASETS, SAVEDATAFILE

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = createdatafile(a,dirr,root,type,par1,par2)

	  
  isdatafile(a);
  if nargin < 6, par2 = []; end
  if nargin < 5, par1 = []; end
  if nargin < 4 || isempty(type), type = 'raw'; end
  if nargin < 3 || isempty(root), root = pwd; end
  if nargin < 2 || isempty(dirr), dirr = getname(a); end
  
  if ~any(strcmp(type,{'raw','cell'}))
    error('Type of datafile is wrong or not supported');
  end
  
  actdir = cd;
  cd(root);
	if exist(dirr) == 7        % if datafile already exists
		cd(dirr);                % go to it
    fclose('all');
		delete *.mat             % make it empty
    delete file*
  else
    [success,mess] = mkdir(dirr);
    if ~success
      error(mess)
    end
    cd(dirr)
  end
    
  if strcmp(type,'raw')
    % construct standard raw datafile
    if isempty(par1), par1 = 'imwrite'; end
    if isempty(par2), par2 = 'bmp'; end
    cmd = par1;
    fmt = par2;
    subdirs = getlablist(a,'string');
    nlab = getnlab(a);
    m = length(nlab);
    s = sprintf('Creating %i files: ',m);
    prwaitbar(m,s);
    n = 0;
    if ~isempty(subdirs)
      for i=1:size(subdirs,1) % run over all classes
        sname = deblank(subdirs(i,:));
        [success,mess] = mkdir(sname);
        if ~success
          error(mess);
        end
        cd(sname)
        J = find(nlab==i);
        p = floor(log10(m))+1;
        form = ['%' int2str(p) '.' int2str(p) 'i'];
        for j = 1:length(J)
          n = n+1;
          prwaitbar(m,n,[s int2str(n)]);
          im = data2im(a(J(j),:));
          feval(cmd,im,['file_' num2str(n,form) '.' fmt],fmt);
        end
        cd ..
      end
    end
    J = find(nlab == 0);
    if ~isempty(J)
      p = floor(log10(m))+1;
      form = ['%' int2str(p) '.' int2str(p) 'i'];
      for j = 1:length(J)
        n = n+1;
        prwaitbar(m,n,[s int2str(n)]);
        im = data2im(a(J(j),:));
        feval(cmd,im,['file_' num2str(n,form) '.' fmt],fmt);
      end
    end
    prwaitbar(0);
    cd ..
    b = prdatafile(dirr);
    if ~isempty(a,'prior')
      b = setprior(b,getprior(a));
    end
    id = getident(a);
    b = setident(b,id);
    cd(dirr);
    save([dirr '.mat'],'b');
    
  elseif strcmp(type,'cell')
    %store datafiles in cells
		if isempty(par1)                   % maximum nr of elements in datafile
			par1 = 1000000;
		end
		if isempty(par2)                   % default: no conversion
			par2 = 0;
		end
    filesize = par1;
    nbits = par2;
  
	                      % reorder datafile for faster retrieval
    file_index = getident(a.prdataset,'file_index');
	  [dummy,R] = sort(file_index(:,1)+file_index(:,2)/(max(file_index(:,2),[],1)+1));
    a = a(R,:);
		[dummy,S] = sort(R); % to reconstruct order later

    dirname = fullfile(root,dirr);
  
    m = length(R);
	  file_index = zeros(m,2);
                             % initialise structure with file information
% skip this for the time being
% cfiles = struct('name',cell(1),'nbits',[], ...
%       'offsetd',[],'scaled',[],'offsett',[],'scalet',[]);
    nobjects = 1;
    files = {};
    s = sprintf('Storing %i objects: ',m);
    prwaitbar(m,s)
		prwaitbar(m,nobjects,[s int2str(nobjects)]);
    im = data2im(a(1,:));
    imsize = numel(im);
    nfiles = 0;
    while (nobjects <= m)
      nfiles = nfiles+1;
      obj1 = nobjects;
      fsize = imsize;
      imcells = {im};
      sizefit = 1;
      mm = 0;
			nobjects = nobjects+1;
			prwaitbar(m,nobjects,[s int2str(nobjects)]);
      while (nobjects <= m) & sizefit
		    im = data2im(a(nobjects,:)); % convert single object
        imsize = numel(im);
        nobjects = nobjects+1;
        mm = mm+1;
        if (fsize+imsize < filesize) & (mm < m/50)
          % store data in same file as long as size fits
          % and never more than 2% of the objects
          % (this is just to keep better track of progress)
          fsize = fsize+imsize;
          imcells = {imcells{:} im};
          sizefit = 1;
          obj2 = nobjects;
        else
          sizefit = 0;
          obj2 = nobjects-1;
					nobjects = obj2;
        end
			end

			obj2 = min(obj2,m);
%		  [imcells,scaled,offsetd,scalet,offsett,prec] = convert(imcells,nbits); % convert precision
      filej = ['file_' num2str(nfiles,'%4.4i')];
		  fname = fullfile(dirname,filej);
      files = {files{:} filej};
      save(fname,'imcells');
%     disp([filej '  ' int2str(length(imcells))]);
%     save(fname,'imcells','scaled','offsetd','scalet','offsett','prec');
    
%     cfiles(j).name    = filej;      % name of the file
%     cfiles(j).prec    = prec;       % precision
%     cfiles(j).scaled  = scaled;     % scaling for datafield
%     cfiles(j).offsetd = offsetd;    % offset for datafield
%     cfiles(j).scalet  = scalet;     % scaling for target field
%     cfiles(j).offsett = offsett;    % offset for target field
% 	  cfiles(j).sized   = a.featsize; % feature (image) size
% 	  cfiles(j).sizet   = size(a.targets,2); % size target field
%                                     % create file_index
		  file_index(obj1:obj2,1) = repmat(nfiles,obj2-obj1+1,1);
		  file_index(obj1:obj2,2) = [1:obj2-obj1+1]';
		end
                                     % finalise datafile definition
    b = prdatafile;                              
    b.files = files;                 % file information
		b.rootpath = dirname;
	  b.type = 'cell';                 %  cell datafile
	  preproc.preproc = [];
	  preproc.pars = {};
	  b.preproc = preproc;             % no preprocessing
	  b.prdataset = setident(a.prdataset,file_index,'file_index'); % store file_index
	  b.prdataset = setname(b.prdataset,dirr);
		b = b(S,:);                      % reconstruct order of objects
	  save(fullfile(dirname,dirr),'b'); % store the datafile as mat-file for fast retrieval
    prwaitbar(0);
  
  end
  cd(actdir);
  
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
	

      
  
return
    