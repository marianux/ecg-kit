function pa = picant(sig,time)
% This function calculates the first peak of signal sig nearest to the end
sig = flipud(sig);
der = diff(sig);
cero = min(find ( (der(1:end-1).*der(2:end)<=0) ) );
% pa_rel = cero+1;
pa = time - cero;  

