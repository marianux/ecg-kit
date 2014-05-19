function ind = buscamin(x);
% It searches the first minimum of the modulus of x.

x = abs(x);
localmin = (x(2:end-1)<=x(1:end-2))...
         & (x(2:end-1)<=x(3:end));

ind = min(find(localmin));
