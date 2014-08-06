%GRIDSIZE Set gridsize used in the plot commands
%
%     O = GRIDSIZE(N)
%
% INPUT
%     N  New grid size (optional, default: display current gridsize)
% 
% OUTPUT
%     O  New grid size (optional)
%		
% DESCRIPTION
% The initial gridsize is 30, enabling fast plotting of PLOTC and PLOTM.
% This is, however, insufficient to obtain accurate graphs, for which a 
% gridsize of at least 100 and preferably 250 is needed. 
% Default: display or return the current gridsize.
%
% EXAMPLES
% See PREX_CONFMAT
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% PLOTC, PLOTM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: gridsize.m,v 1.3 2006/12/19 12:10:06 duin Exp $

function out = gridsize(n)

		persistent CURRENT_GRIDSIZE;

	% If the global variable was not yet initialised, set it to 30 (default).

	if (isempty(CURRENT_GRIDSIZE))
		prwarning(4,'initialising gridsize to 30');
		CURRENT_GRIDSIZE = 30;			
	end

	if (nargin < 1)
		if (nargout == 0)
			disp(['Gridsize is ' num2str(CURRENT_GRIDSIZE) ]);
		end
	else
		if isstr(n)
			n= str2num(n);
		end
		if isempty(n)
			error('Illegal gridsize')
		end
		CURRENT_GRIDSIZE = n;
	end

	if (nargout > 0), out = CURRENT_GRIDSIZE; end

return
