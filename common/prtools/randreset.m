%RANDRESET  Reset state of random generators
%
%   RANDRESET                   - reset states to 1
%   RANDRESET(STATE)            - reset states to STATE
%   STATE = RANDRESET           - retrieve present state (no reset)
%   OLDSTATE = RANDRESET(STATE) - retrieve present state and reset

function output = randreset(state)

if nargin < 1 & nargout < 1
   state = 1;
end

if nargin > 0 
  if iscell(state)
    % old state format, force switch to old generator
    rand('seed',0);
    % new state format, force switch to new generator
  elseif isstruct(state)
    rng(0,'twister');
  end
end

if verLessThan('matlab','8.1.0') || strcmp(getfield(rng,'Type'),'Legacy')
  % old random number generators
  
  if nargout > 0 
    % retrieve present state
    randstate = cell(1,2);
    randstate{1} = rand('state');
    randstate{2} = randn('state');
    output = randstate;
  end

  if ~(nargout == 1 & nargin == 0)
    % reset
    if iscell(state)
      rand('state',state{1});
      randn('state',state{2});
    else
      rand('state',state);
      randn('state',state);
    end
  end

else
  % new random number generator (rng)
  
  if nargout > 0
    % retrieve present state
    output = rng;
  end
  
  if ~(nargout == 1 & nargin == 0)
    % reset
    rng(state);
  end
  
end
  
  
