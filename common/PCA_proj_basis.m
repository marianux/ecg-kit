function autovec = PCA_proj_basis(arg1, arg2, arg3, arg4 )

%sample 10 minutes from the whole recording.
time_sample = 10 * 60; %seconds
time_heartbeat_window = 2; %seconds around the heartbeat

% 10 min of 1000Hz 12-lead ECG
max_samples = 10 * 60 * 12 * 1000; % samples

if( ischar(arg1) && ischar(arg2) )
    
    recording_name = arg1;
    recording_format = arg2;
    q_filters = arg3;
    
    [~, heasig, ann] = read_ECG(recording_name, recording_format );
    
    if( isfield( ann, 'anntyp'))
        ann = AnnotationFilterConvert(ann, recording_format);
    end

    sampling_rate = heasig.freq;
    QRS_locations = ann.time;
    
    cant_QRS_locations = length(QRS_locations);
    %le quito algunos QRS para evitar los de los extremos, y evitar
    %problemas.
    QRS2skip = 10;
    cant_QRS_locations = cant_QRS_locations - QRS2skip;

    cant_QRS_sample = min(cant_QRS_locations, round(time_sample / time_heartbeat_window));

    QRS_sample_idx = randsample(cant_QRS_locations, cant_QRS_sample);
    %corrijo para reubicar los indices segun los QRS q elimine.
    QRS_sample_idx = sort(QRS_sample_idx) + QRS2skip/2;
    
    halfwin_samples = round( time_heartbeat_window/2*heasig.freq);
    some_win = colvec(blackman(2*halfwin_samples+1));
    
    max_block2read = max_samples / heasig.nsig;
    
    ECG_slices = [];
    
    ii = 1;
    while( ii <= cant_QRS_sample)
        this_QRS_sample_idx = QRS_sample_idx( QRS_locations(QRS_sample_idx) >= QRS_locations(QRS_sample_idx(ii)) & QRS_locations(QRS_sample_idx) < QRS_locations(QRS_sample_idx(ii)) + max_block2read );
        this_QRS_locations = QRS_locations(this_QRS_sample_idx);
        this_iter_cant_QRS = length(this_QRS_locations);
        
        this_ECG_start = max(1, this_QRS_locations(1) - round(10*heasig.freq));
        this_ECG_end = min( heasig.nsamp, this_QRS_locations(end) + round(10*heasig.freq));
        
        ECG = read_ECG(recording_name, recording_format, this_ECG_start, this_ECG_end );
        
        %convert ADC samples (int16) to real units data.
        ECG = ADC2realunits(double(ECG), heasig.adczero, heasig.gain);
        
        this_iter_ECG_size = size(ECG,1);
        this_QRS_locations = this_QRS_locations - this_ECG_start + 1;
        
        ECG = BaselineWanderRemovalSplines( ECG, this_QRS_locations, heasig.freq);
        
        aux_idx = arrayfun(@(a)( max(1, this_QRS_locations(a) - halfwin_samples): ...
                                 min( this_iter_ECG_size, this_QRS_locations(a) + halfwin_samples)) , ...
                           colvec(1:this_iter_cant_QRS), 'UniformOutput', false);
        
        ECG_slices = [ECG_slices; cell2mat(cellfun(@(a)( bsxfun(@times, ECG(a,:), some_win) ), aux_idx, 'UniformOutput', false))];
                       
        ii = ii + this_iter_cant_QRS;
    end
    
    ECG = ECG_slices;
    
    clear ECG_slices
    
else
    
    ECG = arg1;
    QRS_locations = arg2;
    sampling_rate = arg3;
    q_filters = arg4;
    
    halfwin_samples = round( time_heartbeat_window/2*sampling_rate);
    some_win = colvec(blackman(2*halfwin_samples+1));
    
    cant_QRS_locations = length(QRS_locations);

    cant_QRS_sample = round(time_sample / time_heartbeat_window);

    QRS_sample_idx = randsample(cant_QRS_locations, cant_QRS_sample);
    
    this_iter_ECG_size = size(ECG,1);

    ECG = BaselineWanderRemovalSplines( ECG, QRS_locations, sampling_rate);

    aux_idx = arrayfun(@(a)( max(1, QRS_locations(QRS_sample_idx(a)) - halfwin_samples): ...
                             min( this_iter_ECG_size, QRS_locations(QRS_sample_idx(a)) + halfwin_samples)) , ...
                       colvec(1:cant_QRS_sample), 'UniformOutput', false);

    ECG = cell2mat(cellfun(@(a)( bsxfun(@times, ECG(a,:), some_win) ), aux_idx, 'UniformOutput', false));
    
end

% Wavelet transform calculation
wtECG = qs_wt(ECG, 4, sampling_rate, q_filters);

% calculate PCA in wt scale 4 of the ECG
result = mcdcov( wtECG, 'plots', 0);
WTecg_cov = result.cov;
clear result
[autovec autoval] = eig(WTecg_cov); 
autoval = diag(autoval);
[~, autoval_idx] = sort(autoval, 'descend');
autovec = autovec(:,autoval_idx);        
