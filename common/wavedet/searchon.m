function onset = searchon(piconset, sig, K)
%
%  onset = searchon(piconset, sig, K);
%
% searches the onset of a wave using the derivative method (I use it with a wavevelet).  sig's last sample must be piconset.
%
% INPUT:
%
%   piconset: position of the first relevant modulos maximum in the wavelet 
%   sig: wavelet signal (one scale)
%   K: threshold factor
%
% Juan Pablo Martínez Cortés
% Last update: 

if isempty(piconset)|isempty(sig), %#ok<OR2>
   onset = [];
   return;
end

maxderon = abs(sig(end));  % maximum derivative
ind1 = min(find(abs(flipud(sig(1:end-1)))<maxderon/K)); %#ok<MXFND>
ind2 = buscamin(flipud(sig(1:end-1)));

if isempty(ind1)&isempty(ind2), %#ok<AND2>
   onset = piconset - length(sig)+1;
elseif isempty(ind1),
   onset = piconset -ind2;
elseif isempty(ind2),
   onset = piconset -ind1;
else
   onset = piconset - min(ind1,ind2);
end

