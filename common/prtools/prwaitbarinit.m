%PRWAITBARINIT Low level routine to simplify PRWAITBAR init
%
%	[N,S,COUNT] = PRWAITBARINIT(STRING,N)
%
% INPUT
%   STRING  - String with text to be written in every waitbar,
%             e.g. 'Processing %i items: '. This will be parsed
%             by S = SPRINTF(STRING,N);
%   N       - Total number of items to be processed
%
% OUTPUT
%   N       - Resulting N (Input N may be expression)
%   S       - Resulting STRING (see above)
%   COUNT   - Counter initialisation, COUNT = 0
%
% This routine has to be used in combination with PRWAITBARNEXT, e.g.:
% 
%  [n,s,count] = prwaitbarinit('Processing %i items:',size(a,1)*size(a,2));
%  for i=1:size(a,1)
%    for j=1:size(a,2)
%      < process a(i,j) >
%      count = prwaitbarnext(n,s,count);
%    end
%  end
%
% PRWAITBARINIT and PRWAITBARNEXT are written on top of PRWAITBAR. If it
% is possible that loops like the above are not fully processed, it is
% necessary to include as a final call PRWAITBAR(0)

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [n,s,count] = prwaitbarinit(ss,n)

s = sprintf(ss,n);
if n > 1
	prwaitbar(n,s);
	prwaitbar(n,1,[s int2str(1)]);
end
count = 0;

return