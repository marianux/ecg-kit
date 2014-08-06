%ISUNTRAINED Test on untrained mapping
%
%    I = ISUNTRAINED(W)
%    ISUNTRAINED(W)
%
% True if the mapping type of W is 'untrained' (see HELP MAPPINGS).
% If called without an output argument ISUNTRAINED generates
% an error if the mapping type of W is not 'untrained'.

% $Id: isuntrained.m,v 1.1 2009/03/18 16:12:41 duin Exp $

function i = isuntrained(w)
	
  if iscell(w)
		i = strcmp(w{1}.mapping_type,'untrained');
		for j=2:length(w)
			if strcmp(w{j}.mapping_type,'untrained') ~= i
				error('Cell array of classifiers cannot be partially trained')
			end
		end
	else
		i = strcmp(w.mapping_type,'untrained');
	end

	if (nargout == 0) & (i == 0)
		error([newline '---- Untrained mapping expected ----'])
	end
return
