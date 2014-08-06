%MAPM Mapping to execute arbitray Matlab function on dataset
%
%    [B,OUT] = MAPM(A,COMMAND,PAR1,PAR2,....)
%    [B,OUT] = A*MAPM([],COMMAND,PAR1,PAR2,....)
%    [B,OUT] = A*MAPM(COMMAND,PAR1,PAR2,....)
%
% INPUT
%    A               Dataset or datafile or double
%    COMMAND         String with function name
%    {PAR1,...}      Cell array with optional parameters for COMMAND
%
% OUTPUT
%    B               Resulting dataset, datafile or double array
%    OUT             (only for datasets or doubles) possible additional
%                    outputs
%
% DESCRIPTION
% On the data in A, the double array A the following command is executed
%
%    [B,OUT] = COMMAND(DATA,PAR1,PAR2,....)
%
% If A is a dataset then B is converted to a similar dataset. If A is a
% cell array, then results are combined into a cell array as well.
%
% This command differs from the similar command FILTM which operates object
% by object (row by row). The COMMAND in MAPM operates on the entire data
% matrix.
%
% EXAMPLES
% argmin = mapm('min')*out2; 
% % A*argmin returns for every column the row indices of the minimum
% % A'*argmin returns for every row the column indices of its minimum
% 
% nexpm = mapm('uminus')*mapm('exp');
% % A*nexpm returns exp(-A), useful incombination with other commands
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES, IM2OBJ, DATA2IM, FILTM, FILTIM, OUT2
% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function varargout = mapm(varargin)

  if nargin == 0
    error('No command found')
  end
  varargin = shiftargin(varargin,'char');
  if isempty(varargin{1})
    varargout{1} = prmapping(mfilename,'fixed_cell',varargin(2:end));
    varargout{1} = setname(varargout{1},'AnyMap');
    return
  end

  % now we have: mapm(data,command,pars)
  % data might be dataset, datafile, cell array, double, structure array
  pars = cell(1,nargin-2);     % input pars
  [a,command,pars{:}] = deal(varargin{:}); % use nice names
  if numel(pars) > 1 % allow argument list as cell array or not.
   pars = {pars}
  end
  varargout = cell(1,nargout);
    
  if isdatafile(a)               % for datafiles processing is stored
        
    if isempty(getpostproc(a)) & ~ismapping(command)
                                 % as preprocessing (if no postproc defined)
      b = addpreproc(a,command,pars);
    else                         % or as mapping as postprocessing
      v = prmapping(mfilename,'fixed',{command,pars});
      v = setname(v,mapname);
      b = addpostproc(a,v);      % store mapping
    end
    varargout{1} = b;
    return
    
  elseif isdataset(a) % are executed here
    
    try
      [varargout{:}] = my_feval(command,a,pars);
    catch me
      [varargout{:}] = my_feval(command,+a,pars);
      if nargout > 0
        if size(varargout{1},1) == size(a,1)
          varargout{1} = setdata(a,varargout{1});
        else
          varargout{1} = prdataset(varargout{1});
          varargout{1} = setuser(varargout{1},getuser(a));
        end
      end
    end
    
  elseif iscell(a)
    
    b = cell(size(a));
    varargout = cell(numel(a),numel(varargout));
    n = numel(b);
    s = sprintf('Processing %i objects: ',n);
    prwaitbar(n,s);
    for i=1:n
      prwaitbar(n,i,[s int2str(i)]);
      [varargout{i,:}] = my_feval(command,a{i},pars);
    end
    prwaitbar(0);
    
  else % doubles
    
    [varargout{:}] = my_feval(command,a,pars{:});
    
  end
  
return

function varargout = my_feval(command,a,pars)
  if nargin < 3, pars = []; end
  varargout = cell(1,nargout);
  if isempty(pars)
    [varargout{:}] = feval(command,a);
  elseif iscell(pars)
    [varargout{:}] = feval(command,a,pars{:});
  else
    [varargout{:}] = feval(command,a,pars);
  end
return