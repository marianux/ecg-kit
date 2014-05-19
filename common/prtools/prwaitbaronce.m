%PRWAITBARONCE Generate single prwaitbar message
%
%	PRWAITBARONCE(STRING,PAR)
%
% INPUT
%   STRING  - String with text to be written in the waitbar,
%             e.g. '%i x %i eigenvalue decomposition: '. 
%             This will be parsed by S = SPRINTF(STRING,PAR{:});
%   PAR     - scalar or cell array with parameter values
%
% This routine has to be used in combination with PRWAITBAR(0), e.g.:
% 
%  prwaitbaronce('%i x %i eigenvalue decomposition ... ',{n,n})
%  [Q,L] = eig(H);
%  prwaitbar(0)
%
% It makes clear to the user what is happening.

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function prwaitbaronce(ss,par)

	if nargin < 2
		prwaitbar(2,ss);
	elseif ~iscell(par)
		prwaitbar(2,sprintf(ss,par));
	else
		prwaitbar(2,sprintf(ss,par{:}));
	end

return