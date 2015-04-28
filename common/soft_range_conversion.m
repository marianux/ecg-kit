%% (Internal) Convert an input range to an output range with a soft function
%
% Sigmoid function mapping an input range to an output range:
% transition_sharpnes: [0-1] 0:  soft trans; 1: hard transition
% 
%     y = soft_range_conversion( x, input_range, output_range, transition_sharpnes )
% 
% 
% Arguments:
% 
%   + x: data within input_range
% 
%   + input_range, output_range: the tanges to convert [min max]
% 
%   + transition_sharpnes: progress_bar handle
% 
% Output:
% 
%   + y: data converted to output_range. 
% 
% Examples:
% 
%       a = linspace(200, 1000, 1000);
%       plot(a, soft_range_conversion(a, [300 600], [0 10], 0.2) )
% 
% 
% See also 
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Birthdate: 17/12/2010
% Last update: 17/12/2010
% Copyright 2008-2015
% 
function y = soft_range_conversion( x, input_range, output_range, transition_sharpnes )

if( nargin < 4 || isempty(transition_sharpnes) )
    transition_sharpnes = 1;
end

a = 10^(5*transition_sharpnes) / diff(input_range);
b = - a * (input_range(1) + diff(input_range)/2);

y = output_range(1) + diff(output_range) * logit_function(a*x+b);
