%ISTRAINED Test on trained mapping
%
%    I = ISTRAINED(W)
%    ISTRAINED(W)
%
% True if the mapping type of W is 'trained' (see HELP MAPPINGS). If
% called without an output argument ISTRAINED generates an error if the
% mapping type of W is not 'trained'.

% $Id: istrained.m,v 1.1 2009/03/18 16:12:41 duin Exp $

function i = istrained(w)

		if iscell(w)
		i = strcmp(w{1}.mapping_type,'trained');
		for j=2:length(w)
			if strcmp(w{j}.mapping_type,'trained') ~= i
				error('Cell array of classifiers cannot be partially trained')
			end
		end
	else
		i = strcmp(w.mapping_type,'trained');
	end

	if nargout == 0 & i == 0
		error([newline '---- Trained mapping expected ----'])
	end

	return
