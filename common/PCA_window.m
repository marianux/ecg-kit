function [autovec autoval] = PCA_window(x, peaks, win_size)
%
% [y,W,A] = PCA_window(x, peaks, win_size)
% PCA in a window centered in peaks(n) of width win_size
%

if( nargin < 3 )
    win_size = 100;
end

win_size = round(win_size/2);

[nsamp, nsig] = size(x);

if( nsamp < nsig )
    x = x';
    [nsamp, nsig] = size(x);
end

peaks = peaks( (peaks - win_size) > 0  & (peaks + win_size ) <= nsamp );
aux_seq_idx = 1:length(peaks);
T0 = cell2mat(arrayfun(@(a)( peaks(a) - win_size:peaks(a) + win_size ), aux_seq_idx(1:end-1), 'UniformOutput', false));

[autovec autoval] = autovec_calculation_robust(x(T0,:));


