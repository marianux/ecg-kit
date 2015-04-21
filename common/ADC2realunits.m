%% Convert adimentional sample values to real units
% Convert adimentional sample values to real units according to zero and
% gain values.
% 
% Example
% 
%   x = ADC2realunits(x, zero, gain)
% 
% See also ADC2units
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Birthdate: 01/01/2012
% Last update: 18/10/2014
% Copyright 2008-2015% Version: 0.1 beta
% Birthdate: 01/01/2012
% Last update: 18/10/2014
% Copyright 2008-2015

function x = ADC2realunits(x, zero, gain)

[CantSamp CantSig CanScales] = size(x);

x = arrayfun( @(a)( bsxfun( @rdivide , bsxfun(@minus, double(squeeze(x(:,:,a))), rowvec(zero)), rowvec(gain))  ), 1:CanScales, 'UniformOutput', false );
x = cell2mat(reshape(x, 1,1, CanScales));

% if( CanScales > 1)
%     x = cellfun(@(a)( bsxfun( @rdivide , bsxfun(@minus, double(a), rowvec(zero)), rowvec(gain)) ), mat2cell(x, CantSamp, CantSig, ones(1,CanScales) ), 'UniformOutput', false);
%     x = cell2mat(x);
% else
%     x = bsxfun( @rdivide , bsxfun(@minus, double(x), rowvec(zero)), rowvec(gain));
% end

