function pp = picpost(sig,time)
% This function calculates the first peak of signal sig nearest to the beginning
der = diff(sig);
cero = min(find ((der(1:end-1).*der(2:end)<=0)));
% pp_rel = ceros+1;
%cero = cero +  (abs(der(cero))>abs(der(cero+1)))
pp = time + cero;  

