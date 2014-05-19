% C = resize(A,k,n,m)
% resizes the n*m image A with border k in 1D vector shape
% into a 2D image C with size n*m and no border
% $Id: resize.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function C = resize(A,k,n,m)
		if nargin < 4
	[n,m] = size(A);
	n = n - 2*k;
	m = m - 2*k;
end
C = reshape(A,n+2*k,m+2*k);
C = C(k+1:k+n,k+1:k+m);
return

