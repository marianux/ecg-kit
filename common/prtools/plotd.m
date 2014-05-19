%PLOTD Plot classifiers, outdated, use PLOTC instead

% $Id: plotd.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function handle = plotd(varargin)

		global PLOTD_REPLACED_BY_PLOTC

	if isempty(PLOTD_REPLACED_BY_PLOTC)
		disp([newline 'PLOTD has been replaced by PLOTC, please use it'])
		PLOTD_REPLACED_BY_PLOTC = 1;
	end

	if nargout > 0
		handle = plotc(varargin{:});
	else
		plotc(varargin{:});
	end

return
