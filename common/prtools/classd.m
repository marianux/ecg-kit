%CLASSD Return labels of classified dataset, outdated, use LABELD instead

% $Id: classd.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function labels = classd(a,w)

		global CLASSD_REPLACED_BY_LABELD

	if isempty(CLASSD_REPLACED_BY_LABELD)
		disp([newline 'CLASSD has been replaced by LABELD, please use it'])
		CLASSD_REPLACED_BY_LABELD = 1;
	end

	if (nargin == 0)

		labels = labeld;

	elseif (nargin == 1)

		labels = labeld(a);

	elseif (nargin == 2)

		labels = labeld(a,w);

	else
		error ('too many arguments');
	end

return

