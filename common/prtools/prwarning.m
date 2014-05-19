%PRWARNING Show PRTools warning
%
%  PRWARNING(LEVEL,FORMAT,...) 
%
% Shows the message (given as FORMAT and a variable number of arguments),
% if the current PRWARNING level is >= LEVEL. Output is written to standard
% error ouput (FID = 2).
%
%  PRWARNING(LEVEL) - Set the current PRWARNING level
%
% Set the PRWARNING level to LEVEL. The default level is 1.
% The levels currently in use are:
%    0  no warnings
%    1  severe warnings (default)
%    2  warnings
%    3  light warnings
%   10  general messages
%   20  program flow messages
%   30  debugging messages
%
%  PRWARNING OFF  - Same as PRWARNING(0)
%  PRWARNING ON   - Same as PRWARNING(1)
%  PRWARNING      - Same as PRWARNING(1)

% Copyright: D. de Ridder, R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: prwarning.m,v 1.5 2007/11/30 16:29:49 duin Exp $

function lev = prwarning (level, varargin)

	persistent PRWARNINGGLOBAL;

	if (isempty (PRWARNINGGLOBAL))
		PRWARNINGGLOBAL = 1;
	end

	if (nargin == 0) & (nargout == 0) % Set warning level to default
		PRWARNINGGLOBAL = 1;
	elseif (nargin == 1)        % Set warning level
		if isstr(level) & strcmp(level,'off')
			level = 0;
		elseif isstr(level) & strcmp(level,'on')
			level = 1;
    elseif isstr(level) % wrong call, no level supplied
      warning(level);   % print standard warning
      level = PRWARNINGGLOBAL; % do not change PRWARNINGGLOBAL
		end
		PRWARNINGGLOBAL = level;
	elseif nargin > 0
		if (level <= PRWARNINGGLOBAL)
			[st,i] = dbstack;   % Find and display calling function (if any)
			if (length(st) > 1)
				caller = st(2).name;
				[paths,name] = fileparts(caller);
				fprintf (2, 'PR_Warning: %s: ', name);
			end;
			fprintf (2, varargin{:});
			fprintf (2, '\n');
		end
	end

	if nargout > 0
		lev = PRWARNINGGLOBAL;
	end
	return
