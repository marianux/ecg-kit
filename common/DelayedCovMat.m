function A = DelayedCovMat(x, peaks, win_size, nsamp, k_hb)
% Internal function for calculating covariance matrix in the vicinity of
% peaks with k_hb heartbeats of delay


    peaks = peaks( (peaks - win_size) > 0  & (peaks + win_size ) <= nsamp );
    aux_seq_idx = 1:length(peaks);
    T0 = cell2mat(arrayfun(@(a)( peaks(a) - win_size:peaks(a) + win_size ), aux_seq_idx(1:end-k_hb), 'UniformOutput', false));
    T1 = cell2mat(arrayfun(@(a)( peaks(a) - win_size:peaks(a) + win_size ), aux_seq_idx(1+k_hb:end), 'UniformOutput', false));
    
    x_T0_centered = bsxfun(@minus,x(T0,:),mean(x(T0,:)));  % Remove mean
    x_T1_centered = bsxfun(@minus,x(T1,:),mean(x(T1,:)));  % Remove mean

    A = 1/length(T0)*(x_T0_centered'*x_T1_centered);

    A = (A+A')/2;
