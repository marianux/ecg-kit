%PRVERSION PRtools version number
%
%		[VERSION,STR,DATE] = PRVERSION
%
% OUTPUT
%		VERSION   Version number (double)
%		STR				Version number (string)
%		DATE			Version date (string)
%
% DESCRIPTION
% Returns the numerical version number of PRTools VER (e.g. VER = 3.0205) 
% and as a string, e.g. STR = '3.2.5'. In DATE, the version date is returned 
% as a string. 

% $Id: prversion.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function [version,str,date] = prversion(location)

if nargin < 1 % get it locally
	signature = prtver;
	str = signature{1}.Version;
	date = signature{1}.Date;
  version = str2version(str);
	%version   = str2num(str(1)) + (str2num(str(3))*1000 + str2num(str(5))*10)/10000;
	if nargout == 0
		disp([newline '   PRTools version ' str newline])
		clear version
  end
else % location is url
  s = urlread(location);
  s = [s ' '];                % make sure there is a space
  s = strrep(s,char(10),' '); % replace newline by space
  n = strfind(lower(s),'version');
  nspace = strfind(s,' ');
  str = s(nspace(1)+1:nspace(2)-1);
  version = str2version(str);
  date = s(nspace(2)+1:nspace(3)-1);
end
  
function version = str2version(str)
n = strfind(str,'.');
if isempty(n)
  version = str2num(str);
elseif length(n) == 1
  version = str2num(str(1:n-1)) + str2num(str(n+1:end))/100;
else
  version = str2num(str(1:n(1)-1)) + str2num(str(n(1)+1:n(2)-1))/100 + str2num(str(n(2)+1:end))/10000;
end
return;
