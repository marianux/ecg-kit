%PRCURSOR Show object ident.
%
%    PRCURSOR(H)
%
% Enable the datacursor in a scatterplot. This can be used to
% investigate the object identifier by clicking on the object.

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function out = prcursor(h,event_obj)

% This file has to tasks: 
% 1. to setup the datacursor (by enabling it and setting the callback
% function)
% 2. to provide the callback when a user pressed on an object

if ~exist('datacursormode','file')
	error('MATLAB 7.0 or newer is required (datacursormode.m is not available).');
end
if nargin ~=2
	% we are doing the setup:
	if nargin<1
		h = gcf;
	end
	dh = datacursormode(h);
	set(dh,'enable','on','updatefcn',@prcursor);

else
	% we are doing the callback:
	nr = get(event_obj,'dataindex');
	ud = get(get(event_obj,'target'),'UserData');
	if ~isempty(ud) & isfield(ud,'ident')
		nr = ud.ident(nr);
	end
	out = sprintf('obj. %d',nr);
end

return
