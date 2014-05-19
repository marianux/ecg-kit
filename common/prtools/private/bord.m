% C = bord(A,n,m)
% Puts a border of width m (default m=1) around image A
% and gives it value n. If n = NaN: mirror image values.
% $Id: bord.m,v 1.3 2007/05/08 15:06:06 duin Exp $

function C = bord(A,n,m);
		if nargin == 2
  m=1;
   prwarning(4,'border width not supplied, assuming 1');
end
[x,y,z] = size(A);
if m > min(x,y)
	mm = min(x,y);
	C = bord(A,n,mm);
	C = bord(C,n,m-mm);
	return
end
if isnan(n)
   C = [A(:,m:-1:1,:),A,A(:,y:-1:y-m+1,:)];
   C = [C(m:-1:1,:,:);C;C(x:-1:x-m+1,:,:)];
else
   bx = ones(x,m,z)*n;
   by = ones(m,y+2*m,z)*n;
   C = [by;[bx,A,bx];by];
end
return
