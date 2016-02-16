%STRLAB Make sure labels are strings
%
%	LAB_OUT = STRLAB(LAB_IN)
%
% The label list LAB_IN (cells numbers or strings) is converted to a
% string matrix LAB_OUT

% $Id: strlab.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function lab2 = strlab(lab1)

if isstr(lab1)
	lab2 = lab1;
elseif iscell(lab1)
	lab2 = char(lab1);
elseif isa(lab1,'double')
	lab2 = num2str(lab1);
else
	error('Illegal label type')
end
