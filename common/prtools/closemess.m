%CLOSEMESS Close progress message
%
%   CLOSEMESS(FID,N)
%
% Closes a progress message of length N on file-id FID
%
% This routine is obsolete now and just preserved to get
% old code running.

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function closemess(fid,n)

	prprogress(fid,'\n');

return
