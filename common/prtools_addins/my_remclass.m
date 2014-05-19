%REMCLASS Remove small classes
%
%   B = REMCLASS(A,N)
%
% INPUT
%   A   Dataset
%   N   Integer, maximum class size to be removed (optional; default 0)
%
% OUTPUT
%   B   Dataset
%
% DESCRIPTION
% Classes having N objects or less are removed. The corresponding objects
% are made unlabeled. Use SELDAT to remove unlabeled objects.
%
% SEE ALSO
% SELDAT, GENDAT

function [b, L] = my_remclass(a,n)
	prtrace(mfilename);
	if nargin < 2, n = 0; end
	
	N = classsizes(a);
	J = find(N <= n);
	L = findnlab(a,J);
	b = setnlab(a,0,L);
	b = setlablist(b);
	
return