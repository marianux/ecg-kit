%NEWLINE The platform dependent newline character
%
% c = newline

% $Id: newline.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function c = newline
	
	if strcmp(computer,'MAC2')
		c = setstr(13);
	elseif strcmp(computer,'PCWIN')
		c = setstr(10);
	else
		c = setstr(10);
	end
	
return
