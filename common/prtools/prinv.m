%PRINV Call to INV() including PRWAITBAR 
%
%	B = PRINV(A)
%
% This calls B = INV(A) and includes a message to PRWAITBAR
% in case of a large A

function B = prinv(A)

[m,n] = size(A);
if min([m,n]) >= 500
	prwaitbaronce('Inverting %i x %i matrix ...',[m,n]);
	B = inv(A);
	prwaitbar(0);
else
	B = inv(A);
end