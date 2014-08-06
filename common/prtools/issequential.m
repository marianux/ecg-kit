%ISSEQUENTIAL Test on sequential mapping
%
%  N = ISSEQUENTIAL(W)
%      ISSEQUENTIAL(W)
%
% INPUT
%  W    input mapping
%
% OUTPUT
%  N    logical value
%
% DESCRIPTION
% Returns true for sequential mappings. If no output is required,
% false outputs are turned into errors. This may be used for
% assertion.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% ISMAPPING, ISSEQUENTIAL

function n = issequential(w)

		
	if isa(w,'prmapping') & strcmp(w.mapping_file,'sequential')
		n = 1;
	else
		n = 0;
	end

	% generate error if input is not a sequential mapping
	% AND no output is requested (assertion)

	if nargout == 0 & n == 0
		error([newline '---- Sequential mapping expected -----'])
	end

return
