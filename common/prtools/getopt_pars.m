%GETOPT_PARS Get optimal parameters from REGOPTC
%
%	   PARS = GETOPT_PARS
%    GETOPT_PARS
%
% DESCRIPTION
% This routine retrieves the parameters as used in the final call
% in computing a classifier they are optimised by REGOPTC.

function pars = getopt_pars

global REGOPT_PARS

if nargout == 0
	s = [];
	for j = 1:length(REGOPT_PARS)
		if isempty(REGOPT_PARS{j})
			s = [s '  []'];
		elseif isstr(REGOPT_PARS{j})
			s = [s '  ' REGOPT_PARS{j}];
		elseif isdataset(REGOPT_PARS{j})
			s = [s '  dataset'];
		elseif ismapping(REGOPT_PARS{j})
			s = [s '  mapping'];
		elseif isstruct(REGOPT_PARS{j})
			s = [s '  structure'];
		else
			s = [s '  ' num2str(REGOPT_PARS{j})];
		end
	end
	fprintf(1,'\nParameters used: %s \n\n',s);
else
	pars = REGOPT_PARS;
end
