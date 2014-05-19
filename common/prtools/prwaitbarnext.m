%PRWAITBARNEXT Low level routine to simplify PRWAITBAR next calls
%
%	COUNT = PRWAITBARNEXT(N,STRING,COUNT)
%
% Update call for PRWAITBAR after initialisation by PRWAITBARINIT.
%
% In case COUNT = 0 it initializes as well and a separate call to
% PRWAITBARINIT is not needed.
%
% SEE PRWAITBARINIT

function count = prwaitbarnext(n,ss,count)

if count < 0
	[n,s,count] = prwaitbarinit(ss,n);
else
	s = sprintf(ss,n);
end
if n > 1
	count = count + 1;
	% we show one count more as prwaitbarnext typically is at the end of a
	% loop, while prwaitbar is at the beginning
	prwaitbar(n,count+1,[s int2str(count+1)]);
	if count == n
		prwaitbar(0);
	end
end
return