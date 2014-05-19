%MAP Map a dataset, train a mapping or classifier, or combine mappings
%
%	B = MAP(A,W) or B = A*W
%
% Maps a dataset A by a fixed or trained mapping (or classifier) W,
% generating
% a new dataset B. This is done object by object. So B has as many objects
% (rows) as A. The number of features of B is determined by W. All dataset
% fields of A are copied to B, except the feature labels. These are defined
% by the labels stored in W.
%
%	V = MAP(A,W) or B = A*W
%
% If W is an untrained mapping (or classifier), it is trained by the dataset A.
% The resulting trained mapping (or classifier) is stored in V.
%
%	V = MAP(W1,W2) or V = W1*W2
%
% The two mappings W1 and W2 are combined sequentially. See SEQUENTIAL for
% a description. The resulting combination is stored in V.
%
% See also DATASETS, MAPPINGS, SEQUENTIAL

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org


function [d,varargout] = map(a,b)
    
%printdebug(a,b) % remove % for printing debug info 
         
if iscell(a) || iscell(b)
  if ismapping(b) && isfixed_cell(b)
    % mappings of type fixed_cell accept cell arrays, treat them as 'fixed'
    b = setmapping_type(b,'fixed');
  else
    % return a cell array with the individual mappings
    % in case and b are both cell arrays, they should have the same size
    % and the result is computed element by element
    d = cellmap(a,b);
    return
  end
end
    
if stamp_map > 0 && nargout == 1 
  % if stamp_map is enabled and ifthis  map is already computed, take it.
  % if map is new and needs to be stored, do so.
	d = stamp_map(a,b);
	if ~isempty(d)
    % this map is found, use it
		return
	end
end
	
% allow for multiple outputs
varargout = repmat({[]},[1, max((nargout-1),0)]);

% enable/disable checking of sizes
global CHECK_SIZES; 
if isempty(CHECK_SIZES), CHECK_SIZES = true; end
% should be false to disable size checking

% get dataset / mapping sizes
[ma,ka] = size(a);
[mb,kb] = size(b);

% batch processing, to be enables/disable by setbatch
% global BATCHSETTINGS % needed to get batch settings in CNORMC
batchsettings = cell(1,3);
d = batchmap(a,b,ma,ka,mb,kb);
if ~isempty(d)
  % in case batch processing is enabled and executed, we are done
  return
end
clear d; % needed to avoid an empty return in case of nargout = 0

if all([ma,mb,ka,kb] ~= 0) && ~isempty(a) && ~isempty(b) && ~isscalar(a) && ~isscalar(b) && ka ~= mb && CHECK_SIZES
	error(['Output size of first argument should match input size of second.' ...
		newline 'Checking sizes might be skipped by defining global CHECK_SIZES = false'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Here the real works starts
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isa(a,'mapping') && isa(b,'mapping')
  % mapping * mapping
	
  if isempty(b) % empty mappings are treated as unity mappings
    d = a;
  elseif istrained(a) && isaffine(a) && istrained(b) && isaffine(b)
		% combine affine mappings
		d = affine(a,b);
  else
    [d, varargout{:}] = sequential(a,b);
	end
	
elseif isa(a,'dataset') || isa(a,'datafile') || isa(a,'double') || ...
    isa(a,'uint8') || isa(a,'uint16') || isa(a,'dipimage') || iscell(a)
  % data * mapping
  
	if isa(a,'uint8') || isa(a,'uint16') || isa(a,'dipimage')
    % convert to double where needed
		a = double(a);
	end
	
	if ~isa(b,'mapping')
		error('Second argument should be mapping or classifier')
  end

  % handle special cases
  if isempty(b) % treat empty mappings as unity mappings
    d = a;
    return
  elseif isempty(a) % return empty mapping if data is empty
		d = mapping([]);
    return
  elseif isa(a,'double') && ka == 1 && ma == 1
	  % handle scalar * mapping by .* (see times.m)
		d = a.*b;
		return; 
  end

  % now we have proper data * mapping
	mapp = getmapping_file(b);

	if isuntrained(b)
    % training
    % set batch settings locally for usage during training, e.g. CNORMC
    % as batch processing might be needed inside training
		pars = +b; % data stored in mapping
    % might be another mapping
    if issequential(b) || isstacked(b) || isparallel(b)
      [d, varargout{:}] = feval(mapp,a,b); 
		else
			if ~iscell(pars), pars = {pars}; end
      [d, varargout{:}] = feval(mapp,a,pars{:});
    end
    if ~isa(d,'mapping')
			error('Training an untrained classifier should produce a mapping')
    end
    
    % training is finished, some postprocessing needed
    if getout_conv(b) > 1
      % proper confidences needed, set out_conv in mapping
			d = d*classc;
    end
    % set scaling
		d = setscale(d,getscale(b)*getscale(d));
		name = getname(b);
    if ~isempty(name)
			d = setname(d,name);
    end
    [batchsettings{:}] = getbatch(b);
    d = setbatch(d,batchsettings{:}); % copy batch from untrained classifier
    
	elseif isdatafile(a) && istrained(b)
    % execution of datafile on trained mapping
 		if issequential(b)
 			d = feval(mapp,a,b);
 		else
			d = addpostproc(a,{b}); % just add mapping to postprocesing and execute later
		end
		
	elseif isdatafile(a)
    % execution of datafile on fixed mapping
		try  % try whether this is a mapping that knows how to handle a datafile
			pars = getdata(b); % parameters supplied in fixed mapping definition
			if nargout > 0
        if isempty(varargout)
					d = feval(mapp,a,pars{:});
				else    
					[d, varargout{:}] = feval(mapp,a,pars{:});
        end 
        if b.scale ~= 1
          d = b.scale*d;
        end
			else
				b.scale*feval(mapp,a,pars{:})
				return
			end
    catch
			[lastmsg,lastid] = lasterr;
      if ~strcmp(lastid,'prtools:nodatafile')
        % rethrow is buggy, so generate again the erroneous call
        disp(lastid)
        feval(mapp,a,pars{:});
				error(lastmsg);
      end
      d = addpostproc(a,{b}); % just add mapping to postprocesing and execute later
      return
		end
  elseif isfixed(b) & isparallel(b)
    % execution of fixed paralllel mapping
    d = parallel(a,b);
 	elseif isfixed(b) || isfixed_cell(b) || iscombiner(b)
    % execution of fixed or combiner
		pars = getdata(b); % parameters supplied in fixed mapping definition
    if ~iscell(pars), pars = {pars}; end
		if nargout == 1
			fsize = getsize_in(b);
			if isdataset(a) & any(fsize~=0) & ~isobjim(a)
				a = setfeatsize(a,fsize); % needed to set object images, sometimes
			end 
			d = feval(mapp,a,pars{:});
    elseif nargout > 1    
      [d, varargout{:}] = feval(mapp,a,pars{:});
    else % no output parameters, display result
			feval(mapp,a,pars{:})
			return
		end

	elseif istrained(b)
    % execution of trained mapping
		%if ~isdataset(a)  % needed ?????
		%	a = dataset(a);
		%end
    if isa(a,'dataset')
			fsize = getsize_in(b);
			if any(fsize~=0) & ~isobjim(a)
				a = setfeatsize(a,fsize); % needed to set object images, sometimes
      end
    end
		[d,varargout{:}] = feval(mapp,a,b);
    if isa(a,'dataset')
      d = setuser(d,getname(b),'evaluated_by');
    end
    if ~isreal(+d)
			prwarning(2,'Complex values appeared in dataset');
    end
    if ~isdataset(d) && getout_conv(b) > 0
      % classifiers should output dataset to set featue labels (classnames)
      d = setfeatlab(dataset(d),getlabels(b));
    end
		if isdataset(d)
			d = setcost(d,b.cost);
			% see if we have reasonable data in the dataset
		end

	else
		error(['Unknown mapping type: ' getmapping_type(b)])
	end
	
  if isdataset(d) 
	
			% we assume that just a basic dataset is returned, but that scaling
			% and outputconversion still have to be done.

			% scaling
		v = getscale(b);
		if length(v) > 1, v = repmat(v(:)',ma,1); end
		d = v.*d;
			% outputconversion
		switch 	getout_conv(b);
		case 1  % SIGM output
			if size(d,2) == 1
				d = [d -d]; % obviously still single output discriminant
				d = setfeatlab(d,d.featlab(1:2,:));
			end             
			d = sigm(d);
		case 2  % NORMM output
			if size(d,2) == 1
				d = [d 1-d]; % obviously still single output discriminant
				d = setfeatlab(d,d.featlab(1:2,:));
			end             % needs conversion to two-classes before normm
			d = normm(d);
		case 3  % SIGM and NORMM output
			if size(d,2) == 1
				d = [d -d]; % obviously still single output discriminant
				d = setfeatlab(d,d.featlab(1:2,:));
			end             % needs conversion to two-classes before sigmoid
			d = sigm(d);
			d = normm(d);
		end
		%DXD finally, apply the cost matrix when it is defined in the
		%mapping/dataset:
		d = costm(d);
	end

elseif isa(a,'mapping')
  % mapping * data
	if isa(b,'dataset')
		error('Datasets should be given as first argument')
	elseif isdouble(b) && isscalar(b)
		d = setscale(a,b*getscale(a));
	elseif istrained(a) && isdouble(b)
		d = a*affine(b);
	else
		error('Mapping not supported')
	end
		
else
	%a
	b
	error('Data type not supported')
end

function d = batchmap(a,b,ma,ka,mb,kb)

% compute d in batch mode, return d = [] if not applicable

d = [];
if ismapping(b) && ~isuntrained(b) % check for batch processing
  [batchflag,batchsize,objsize] = getbatch(b);
  if batchflag && ((isdataset(a) && ~isfeatim(a)) || isdouble(a))  && (ma > objsize)
    % go for batch processing
    s = sprintf('Mapping %i objects: ',ma);
    prwaitbar(ma,s);
    %DXD map the first batch to setup the output dataset:
    dd = map(a(1:batchsize,:),b);
    %DXD first test if we are dealing with a mapping that outputs just a
    %single value (like testc):
    nb = size(dd,1);
    average_output = 0;
    if (nb~=batchsize)
      if (nb==1) % map returns a single value
        warning('prtools:map:AverageBatchOutputs',...
        ['The mapping appears to return a single object from a input',...
        newline,...
        'dataset. The objects resulting from different batches in the batch',...
        newline,'processing will be *averaged*.']);
        average_output = 1;
      end
    end
    kb = size(dd,2);
    nobatch = false;
    if isdataset(dd)
      if isdataset(a)
        d = setdata(a,zeros(ma,kb));
      else
        d = dataset(zeros(ma,kb));
      end
      d = setfeatlab(d,getfeatlab(dd));
    elseif isa(dd,'double')
      d = zeros(ma,kb);
    else 
		% irregular, escape from batch processing
      d = [];
      prwaitbar(0);
      return
    end
    d(1:batchsize,:) = dd;
    
      % start the batch processing
      n = floor(ma/batchsize);
      prwaitbar(ma,batchsize,[s int2str(batchsize)]);
      for j=2:n
        L = (j-1)*batchsize+1:j*batchsize;
        aa = doublem(a(L,:));
        d(L,:) = map(aa,b);
        prwaitbar(ma,j*batchsize,[s int2str(j*batchsize)]);
      end
      L = n*batchsize+1:ma;
      if ~isempty(L)
        aa = doublem(a(L,:));
        dd = map(aa,b);
        d(L,:) = dd;
      end
      if isdataset(d)
        featlabd = getfeatlab(d);
        d = setdat(a,d);
        d = setfeatlab(d,featlabd);
      end
      if average_output
        d = mean(d,1);
      end
      prwaitbar(0);
  end
end


function d = cellmap(a,b)

if iscell(a) && ~iscell(b)
  d = cell(size(a));
  [n,s,count] = prwaitbarinit('Mapping %i cells: ',numel(a));
  for i = 1:size(a,1);
  for j = 1:size(a,2);
    d{i,j} = map(a{i,j},b);
    count = prwaitbarnext(n,s,count);
  end
  end
elseif ~iscell(a) && iscell(b)
  d = cell(size(b));
  [n,s,count] = prwaitbarinit('Mapping %i cells: ',numel(b));
  for i = 1:size(b,1);
  for j = 1:size(b,2);
    d{i,j} = map(a,b{i,j});
    count = prwaitbarnext(n,s,count);
  end
  end
else
  d = cell(numel(a),numel(b));
  [n,s,count] = prwaitbarinit('Mapping %i cells: ',numel(d));
  for i = 1:length(a)
    for j = 1:length(b)
      d{i,j} = map(a{i},b{j});
      count = prwaitbarnext(n,s,count);
    end
  end
end

function printdebug(a,b)
% print debug info
  disp(' ')
  if isdouble(a)
    disp(['double ' num2str(size(a,1)) ' by ' num2str(size(a,2))])
  else
    a
  end
  if isdouble(b)
    disp(['double ' num2str(size(b,1)) ' by ' num2str(size(b,2))])
  else
    b
  end
return
    
  

