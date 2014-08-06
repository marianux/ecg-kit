%PROCM Fixed-cell mapping to execute arbitray Matlab function on objects
%      Low-level routine
%
%    [OUT1,OUT2, ..] = A*PROCM(COMMAND,PAR1,PAR2,....)
%    [OUT1,OUT2, ..] = {A,B}*PROCM(COMMAND,PAR1,PAR2,....)
%
% INPUT
%    A               Dataset, datafile, cell array or double
%    B               Dataset, datafile, cell array or double
%    COMMAND         String with function name or PRTools nmapping
%    PAR1,...        Optional parameters for COMMAND
%
% OUTPUT
%    OUT1, ...       Cell arrays with results from COMMAND
%
% DESCRIPTION
% This routine executes a given COMMAND over a set of objects. These can
% be supplied by a dataset, a datafile, a cell array or a double array.
% Rows are interpreted as objects. For each object N stored in A the
% routine COMMAND is processed by
%
%    [OUT1{N}, ...] = COMMAND(OBJECT_N,PAR1,PAR2,....)
%
% For diadic operations the two inputs should be combined in a cell: {A,B}.
% They are processed by
%
%    [OUT1{N}, ...] = COMMAND(OBJECT_A_N,OBJECT_B_N,PAR1,PAR2,....)
%
% This is a low-level routine that feeds the unpacked objects in datafiles
% and dataset (doubles) to COMMAND. Use FILTM for a higher level of
% processing. Use MAPM to handle entire matrices, datasets or datafiles
% instead of processing object by object.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES, IM2OBJ, FEAT2OBJ, DATA2IM, FILTM, FILTM, MAPM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function varargout = procm(varargin)

  if nargin == 0
    error('No command found')
  end
  varargin = shiftargin(varargin,'char');
  if isempty(varargin{1})
    w = prmapping(mfilename,'fixed_cell',varargin(2:end));
    varargout{1} = setname(w,'DataProc');
    return
  end

  % now we have: procm(data,command,pars)
  % data might be dataset, datafile, cell array, double, structure array
  
  nout = max(1,nargout);
  varargout = cell(1,nout); % output list
  argout = cell(1,nout);    % temp cell array for output per object
  pars = cell(1,nargin-2);     % input pars
  [a,command,pars{:}] = deal(varargin{:}); % use nice names
%   if numel(pars) > 1 % allow argument list as cell array or not.
%    pars = {pars}
%   end
  if iscell(a) && all(size(a)==[1,2])      % diadic operation
    [a,b] = deal(a{:});
    diadic = true;
  else
    diadic = false;
  end
  m = size(a,1);
  for j=1:nout              % output space
    varargout{j} = cell(m,1);
  end
  t = sprintf(['Processing %i objects by ' command ': '],m);
  prwaitbar(m,t)
  for i=1:m                   % loop over objects
    prwaitbar(m,i,[t num2str(i)]);
    if diadic                 % call command
      [argout{:}] = feval(command,g(a,i),g(b,i),pars{:});
    else
      x = g(a,i);
      [argout{:}] = feval(command,g(a,i),pars{:});
    end
    for j=1:nout           % store results
      varargout{j}{i} = argout{j}; 
    end
  end
  prwaitbar(0);

 
return

function x = g(a,i)
  if iscell(a)              % get cell data
    x =a{i}; 
  elseif isa(a,'prdataset') % datasets and datafiles
    x = +a(i,:);            % get their data
  else                      % all other datatypes
    x = a(i,:);
  end
return
