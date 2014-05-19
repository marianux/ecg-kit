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

if nargout > 0
	randstate = cell(1,2);
	randstate{1} = rand('state');
	randstate{2} = randn('state');
	output = randstate;
end

if ~(nargout == 1 & nargin == 0)
	if iscell(state)
		rand('state',state{1});
		randn('state',state{2});
	else
		rand('state',state);
		randn('state',state);
	end
end

   
