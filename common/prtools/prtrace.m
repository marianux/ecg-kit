%PRTRACE Trace PRTools routines
%
% Routine is outdated and will directly return
%
%	PRTRACE ON              Tracing of the PRTools routines is switched on
% PRTRACE OFF             Tracing of the PRTools routines is switched off
%	PRTRACE(MESSAGE,LEVEL)  
%
% INPUT
%   MESSAGE  String 
%   LEVEL    Trace level (optional; default: 1/0 if PRTRACE is ON/OFF)   
%
% DESCRIPTION
% Note that PRTRACE describes both a global variable (ON/OFF) and a function.
% If PRTRACE is ON, each access to each PRTools routine is reported. 
% Additional trace messages may be added by PRTRACE(MESSAGE) in which MESSAGE
% is a string, listed when the command is encountered and PRTRACE is ON. The 
% standard MESSAGE is the name of the executing m-file (MFILENAME).
%
%	PRTRACE(MESSAGE,LEVEL) lists the traced message above some level. If 
% PRTRACE is switched ON, the trace-level is increased by one. Only 
% if LEVEL >= trace-level, MESSAGE is displayed. If PRTRACE is switched 
% OFF, the trace-level is reset to 0. For all main routines of PRTools, 
% LEVEL is set to 1. All basic routines defining MAPPINGS and DATASETS 
% (except for MAPPING and DATASET themselves) have the trace-level of 2.

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: prtrace.m,v 1.4 2008/07/28 09:02:30 duin Exp $

function prtrace(message,level)

return
% 	persistent PRTRACEGLOBAL
% 
% 	% If PRTRACEGLOBAL was not defined, do it now.
% 	if isempty(PRTRACEGLOBAL)
% 		PRTRACEGLOBAL = 0;
% 		dataset;   % make sure that PRTools is activated
% 	end
% 
% 	% Check the input arguments.
% 	if (nargin == 0)
% 		;
% 	elseif (~isstr(message))
% 		error('No MESSAGE string found.')
% 	else
% 		% Depending on message, we should turn the messaging ON or OFF, or we
% 		% have to increase the trace-level and show the message when the
% 		% trace-level is high enough.
% 
% 		switch message
% 		case {'on','ON'}
% 			PRTRACEGLOBAL = PRTRACEGLOBAL+1;
% 		case {'off','OFF'}
% 			PRTRACEGLOBAL = 0;
% 		otherwise
% 			if (PRTRACEGLOBAL)
% 				if (nargin == 1) 
% 					level = 1;
% 				end
% 		 		if (PRTRACEGLOBAL >= level)
% 					disp(message);
% 				end
% 			end
% 		end
% 	end
% 
% return;
