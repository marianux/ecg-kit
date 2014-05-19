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

