%% (Internal) Pattern matching function to be used in an arbitrary task
%
%     similarity = similarity_calculation(signal, ECG_header_aux, progress_handle, pattern2detect )
% 
% 
% Arguments:
% 
%   + signal: 
% 
%   + ECG_header_aux: header info
% 
%   + progress_handle: progress_bar handle
% 
%   + pattern2detect: the pattern to detect in signal
% 
% Output:
% 
%   + similarity: a signal proportional to the similarity between
%   pattern2detect and signal. 
% 
% Example:
% 
% 
% See also QRScorrector, ECGtask_arbitrary_function
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Birthdate: 17/12/2010
% Last update: 17/12/2010
% Copyright 2008-2015
% 
function similarity = similarity_calculation(signal, ECG_header_aux, progress_handle, pattern2detect )

    [nsamples_signal, nsig_signal] = size(signal);
    [nsamples_pattern2detect, nsig_pattern2detect] = size(pattern2detect);
    
    if( nsig_signal ~= nsig_pattern2detect )
        error('similarity_calculation:BadArg', 'Signal and pattern to seek must have the same amount of columns.')
    end
    
    similarity = cellfun( @(a,b)( diff(conv( a, b, 'same' )) ), mat2cell(double(signal), nsamples_signal, ones(1,nsig_signal)), mat2cell(pattern2detect, nsamples_pattern2detect, ones(1,nsig_signal)), 'UniformOutput', false);
    similarity = cell2mat(cellfun( @(a,b)( diff(conv( [a(1,:);a], b, 'same' )) ), similarity, mat2cell(flipud(pattern2detect), nsamples_pattern2detect, ones(1,nsig_signal)), 'UniformOutput', false));
    similarity = (mean( [similarity(1,:);similarity] ,2));
    