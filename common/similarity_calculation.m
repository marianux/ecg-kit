function similarity = similarity_calculation(signal, ECG_header_aux, progress_handle, pattern2detect )

    [nsamples_signal, nsig_signal] = size(signal);
    [nsamples_pattern2detect, nsig_pattern2detect] = size(pattern2detect);
    
    if( nsig_signal ~= nsig_pattern2detect )
        error('similarity_calculation:BadArg', 'Signal and pattern to seek must have the same amount of columns.')
    end
    
    similarity = cellfun( @(a,b)( diff(conv( a, b, 'same' )) ), mat2cell(double(signal), nsamples_signal, ones(1,nsig_signal)), mat2cell(pattern2detect, nsamples_pattern2detect, ones(1,nsig_signal)), 'UniformOutput', false);
    similarity = cell2mat(cellfun( @(a,b)( diff(conv( [a(1,:);a], b, 'same' )) ), similarity, mat2cell(flipud(pattern2detect), nsamples_pattern2detect, ones(1,nsig_signal)), 'UniformOutput', false));
    similarity = (mean( [similarity(1,:);similarity] ,2));
    