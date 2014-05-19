%RUNCELLS Run command for all cell elements of first parameter

function varargout = runcells(command,varargin)

first = varargin{1};
argout = cell(1,nargout);
seed = randreset;
varargout = cell(1,nargout);
if ~iscell(first)
  [varargout{:}] = feval(command,varargin{:});
else
  for i=1:nargout
    varargout{i} = cell(size(first));
  end

  for j=1:numel(first)
    randreset(seed)
    if nargin == 2
      [argout{:}] = feval(command,first{1});
    else
      [argout{:}] = feval(command,first{1},varargin{2:end});
    end
    for i=1:nargout
      varargout{i}{j} = argout{i};
    end
  end
end
return

