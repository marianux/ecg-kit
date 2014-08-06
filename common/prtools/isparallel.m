%ISPARALLEL Test on parallel mapping
%
%  N = ISPARALLEL(W)
%      ISPARALLEL(W)
%
% INPUT
%  W    input mapping
%
% OUTPUT
%  N    logical value
%
% DESCRIPTION
% Returns true for parallel mappings. If no output is required,
% false outputs are turned into errors. This may be used for
% assertion.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% ISMAPPING, ISSTACKED

% $Id: isparallel.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function n = isparallel(w)

		
	if isa(w,'prmapping') & strcmp(w.mapping_file,'parallel')
		n = 1;
	else
		n = 0;
	end

	% generate error if input is not a parallel mapping
	% AND no output is requested (assertion)

	if nargout == 0 & n == 0
		error([newline '---- Parallel mapping expected -----'])
	end

return
