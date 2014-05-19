function offset = searchoff(picoffset, sig, K)
%
% offset = searchoff(picoffset, sig, K);
%
% searches the offset of a wave using the derivative method (I use it with a wavevelet).  sig's first sample must be piconset.
%
% INPUT:
%
%   picoffset: position of the last relevant modulos maximum in the wavelet 
%   sig: wavelet signal (one scale)
%   K: threshold factor
%
% Juan Pablo Martínez Cortés
% Last update: 

if isempty(picoffset)|isempty(sig), %#ok<OR2>
   offset = [];
   return;
end

maxderoff = abs(sig(1));  % maximum derivative
ind1 = min(find(abs(sig(2:end))<maxderoff/K)); %#ok<MXFND>
ind2 = buscamin(sig(2:end));

if isempty(ind1)&isempty(ind2), %#ok<AND2>
   offset = picoffset + length(sig)-1;
elseif isempty(ind1),
   offset = picoffset +ind2;
elseif isempty(ind2),
   offset = picoffset +ind1;
else
   offset = picoffset +min(ind1,ind2);
end

