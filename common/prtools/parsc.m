%PARSC Parse classifier
% 
% 	PARSC(W)
% 
% Displays the type and, for combining classifiers, the structure of the
% mapping W.
% 
% See also MAPPINGS

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: parsc.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function parsc(w,space)

		% If W is not a mapping, do not generate an error but simply return.

	if (~isa(w,'prmapping')) return; end
	if (nargin == 1)
		space = ''; 
	end
	
	% Display W's type.

	display(w,space);

	% If W is a combining classifier, recursively call PARSC to plot the
	% combined classier.

	pars = getdata(w);
	if (iscell(pars))
		space = [space '  '];
		for i=1:length(pars)
			parsc(pars{i},space)
		end
	end

return

