%PRPROGRESS Report progress of some PRTools iterative routines
%
%  PRPROGRESS ON
%
% All progress of all routines will be written to the command window.
%
%  PRPROGRESS(FID)
%
% Progress reports will be written to the file with file descriptor FID.
%
%  PRPROGRESS(FID,FORMAT,...)
%
% Writes progress message to FID. If FID == [], the predefined destination
% (command window or file) is used. 
%
%  PRPROGRESS OFF
%  PRPROGRESS(0)
%
% Progress reporting is turned off.
%
%  PRPROGRESS
%
% Toggles between PRPROGRESS ON and PRPROGRESS OFF
%
%  FID = PRPROGRESS
%
% Retrieves the status of PRPROGRESS
%
% Some routines (e.g. CLEVAL)  have a switch in the function call by which
% progress reporting for that routine only can be initiated.
%
% By default, PRPROGRESS is switched off. Interactive progress tracing can
% be best following by PRWAITBAR
%
% SEE ALSO
% PRWAITBAR

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function fid = prprogress(par,varargin)

	persistent GLOBALPRPROGRESS
  mlock; % prevents GLOBALPRPROGRESS being cleared by CLEAR ALL

	if isempty(GLOBALPRPROGRESS)
		GLOBALPRPROGRESS = 0;
	end

	if nargin == 0 
		
		if nargout == 0
			if GLOBALPRPROGRESS ~= 0, 
				GLOBALPRPROGRESS = 0;
			else
				GLOBALPRPROGRESS = 1; 
			end
		else
			fid = GLOBALPRPROGRESS;
		end

	elseif nargin == 1 & nargout == 0
		
		if isstr(par)
			switch par
			case {'on','ON'}
				GLOBALPRPROGRESS = 1;
			case {'off','OFF'}
				GLOBALPRPROGRESS = 0;
			otherwise
				error('Illegal input for PRPROGRESS')
			end
		else
			GLOBALPRPROGRESS = par;
		end
		
	elseif (GLOBALPRPROGRESS > 0) & (~isempty(par)) & (par>0)  %DXD, fid=0 means no printing...
		
		s = sprintf(varargin{:});
		n = fprintf(par,s);
		
		if nargout > 0
			fid = length(s);
		end
		
	elseif GLOBALPRPROGRESS > 0
		
		s = sprintf(varargin{:});
		n = fprintf(GLOBALPRPROGRESS,s);
		
		if nargout > 0
			fid = length(s);
		end
		
	else
		
		if nargout > 0
			fid = 0;
		end
		
	end
	
return
	

