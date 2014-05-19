%PRPINV Call to PINV() including PRWAITBAR 
%
%	B = PRPINV(A)
%
% This calls B = PINV(A) and includes a message to PRWAITBAR
% in case of a large A

function B = prpinv(A)

[m,n] = size(A);
if min([m,n]) > 250
	prwaitbaronce('Inverting %i x %i matrix ...',[m,n]);
	B = pinv(A);
	prwaitbar(0);
else
	B = pinv(A);
end