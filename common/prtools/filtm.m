%FILTM Mapping to filter objects in datasets and datafiles
%
%    B = FILTM(A,COMMAND,{PAR1,PAR2,....},SIZE)
%    B = A*FILTM([],COMMAND,{PAR1,PAR2,....},SIZE)
%    B = A*FILTM(COMMAND,{PAR1,PAR2,....},SIZE)
%
% INPUT
%    A               Dataset or datafile or double
%    COMMAND         String with function name of command to be executed
%    {PAR1,...}      Cell array with optional parameters to COMMAND
%    SIZE            Output size of the mapping (default: input size)
%
% OUTPUT
%    B               Dataset or datafile of images processed by COMMAND
%
% DESCRIPTION
% For each object stored in A a filter operation is performed as
%
%    OBJECT_OUT = COMMAND(OBJECT_IN,PAR1,PAR2,....)
%
% The results are collected and stored in B. In case A (and thereby B) is
% a datafile, execution is postponed until conversion into a dataset, or a
% call to SAVEDATAFILE or CREATEDATAFILE.
%
% EXAMPLES
% B = A*filtm('histc',[1:256]);
% Computes for every object in A a histogram and stores it in B.
%
% A similar command FILTIM is recommended for handling multiband images.
% The more general command PROCM is useful for obtaining multiple outputs.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES, IM2OBJ, DATA2IM, IM2FEAT, DATGAUSS, DATFILT, FILTIM
% SAVEDATAFILE, CREATEDATAFILE, MAPM, PROCM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function b = filtm(varargin)

  if nargin == 0
    error('No command found')
  end
  
  % check for mapping definition: first argument empty or string
  varargin = shiftargin(varargin,'char');
  if isempty(varargin{1})
    b = prmapping(mfilename,'fixed_cell',varargin(2:end));
    b = setname(b,'DataFilt');
    % set size_out if user wants so
    if numel(varargin) == 4, b = setsize_out(b,varargin{4}); end
    return
  end

  % now we have: proc(data,command,pars)
  % data might be dataset, datafile, cell array, double, structure array
  nout = max(1,nargout);
  argin = setdefaults(varargin,[],[],[],[]);
  [a,command,pars,outsize] = deal(argin{:}); % use nice names
  if isempty(pars)
    pars = {};
  elseif ~iscell(pars)
    pars = {pars}; 
  end
    
  if isdatafile(a)               % for datafiles filters are stored
        
    if isempty(getpostproc(a)) && ~ismapping(command)
                                 % as preprocessing (if no postproc defined)
      b = addpreproc(a,command,pars,outsize);
    else                         % or as mapping as postprocessing
      if ismapping(command)      % we have already a mapping
        v = command;
      else                       % just a string, construct mapping
        v = prmapping(mfilename,'fixed',{command,pars});
        v = setname(v,mapname);
      end
      if ~isempty(outsize)       % user wants to set an output size (hope he has good reasons)
        v = setsize_out(v,outsize); % add it to mapping
      end
      b = addpostproc(a,v);      % store mapping
    end
    if ~isempty(outsize)
      b = setfeatsize(b,outsize); % set featsize equal to output size
    end
    return
    
  elseif isdataset(a) % datasets are executed here
   
    m = size(a,1);                 
    d = +a;
    imsize = getfeatsize(a);    
    if length(imsize) == 1
      imsize = [1 imsize];
    end
    
    % Perform command on first object to check whether image size stays equal
    first = execute(command,reshape(d(1,:),imsize),pars);
    first = double(first); % for DipLib users
    if isempty(outsize)
      outsize = size(first);
    end
    
    % process all other objects
    out = repmat(first(:)',m,1);
    for i = 2:m
      ima = double(execute(command,reshape(d(i,:),imsize),pars));
      sima = size(ima);
      if (any(outsize ~= sima(1:length(outsize))))
        error('All object sizes should be the same')
      end
      out(i,:) = ima(:)';
    end
    
    % store processed objects in dataset
    b = setdata(a,out);
    b = setfeatsize(b,outsize);
    
  elseif iscell(a)
    
    b = cell(size(a));
    n = numel(b);
    s = sprintf('Filtering %i objects: ',n);
    prwaitbar(n,s);
    for i=1:n
      prwaitbar(n,i,[s int2str(i)]);
      b{i} = feval(mfilename,a{i},command,pars,outsize);
    end
    prwaitbar(0);
    
  else
    
    b = +feval(mfilename,prdataset(a),command,pars);
    
  end
  
return

function out = execute(command,a,pars)

  exist_command = exist(command);
  if isstr(command) & any([2 3 5 6] == exist_command)
    if isempty(pars)
      out = feval(command,a);
    else
      out = feval(command,a,pars{:});
    end
  elseif ismapping(command)
    out = prmap(a,command);
  else
    error('Filter command not found')
  end
 
return