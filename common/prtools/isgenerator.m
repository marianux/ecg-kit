%ISGENERATOR Test on generator mapping
%
%    I = ISGENERATOR(W)
%    ISGENERATOR(W)
%
% True if the mapping type of W is 'generator' (see HELP MAPPINGS). If 
% called without an output argument ISGENERATOR generates an error if the 
% mapping type of W is not 'generator'.


function i = isgenerator(w)

		
	i = strcmp(w.mapping_type,'generator');

	if nargout == 0 && i == 0
		error([newline '---- Generator mapping expected ----'])
	end

	return
