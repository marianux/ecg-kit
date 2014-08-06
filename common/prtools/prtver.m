%PRTVER Get PRTools version
%
%This routine is intended for internal use in PRTools only

function prtversion = prtver

persistent PRTVERSION
if ~isempty (PRTVERSION)
	prtversion = PRTVERSION;
else
  [dummy,prtoolsname] = fileparts(fileparts(which('fisherc')));
  prtversion = {ver(prtoolsname) datestr(now)};
  PRTVERSION = prtversion;
end
