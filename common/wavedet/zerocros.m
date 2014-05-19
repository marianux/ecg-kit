function index = zerocros(x),
% Returns the index of the input vector in which the first zero crossing is located.
% Juan Pablo Martínez Cortés
% Last update:  
%
% tested with MATLAB Version R13
%

m = x(2:end).*x(1:end-1);
index = min(find(m<=0));
if abs(x(index))>abs(x(index+1)),
  index = index+1;
end



