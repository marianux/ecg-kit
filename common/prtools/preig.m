%PREIG Call to EIG() including PRWAITBAR 
%
%	[E,D] = PREIG(A)
%
% This calls [E,D] = EIG(A) and includes a message to PRWAITBAR
% in case of a large A

function [E,D] = preig(A)

[m,n] = size(A);
if min([m,n]) > 500
	prwaitbaronce('Computing %i x %i eigenvectors ...',[m,n]);
	if nargout < 2
		E = eig(A);
	else
		[E,D] = eig(A);
	end
	prwaitbar(0);
else
	if nargout == 1
		E = eig(A);
	else
		[E,D] = eig(A);
	end
end