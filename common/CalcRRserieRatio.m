function [ ratio, estimated_labs ]= CalcRRserieRatio(time_serie, ECG_header, start_end)

% Description: 
% 
% Calculates the quality of the QRS complex detections in time_serie. 
% 
% Arguments:
%     
%     +time_serie: [cell] REQUIRED
%           
%           Cell array of size [nsig 1] with the detections
%           performed for each lead.
% 
%     +ECG_header: [struct] OPTIONAL. 
% 
%             Description of the ECG typically available in the
%             ECG_header. Structure with fields:
% 
%               -freq: Sampling rate in Hz. (1)
% 
%               -nsig: Number of ECG leads. (size(ECG,2))
% 
%               -nsamp: Number of ECG samples. (size(ECG,1))
%                 
%     +start_end: [numeric] OPTIONAL. Vector of size [2 1] with the start
%                           and end indexes.
%     
% 
% Limits and Known bugs:
%   Probably a lot :( ... but dont panic! send me feedback if you need help.
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 23/4/2013

    
    %% constants

    min_pattern_separation = 350;  % ms
    max_pattern_separation = 1500; % ms

    %% start

    if( nargin < 2 || isempty(time_serie) || isempty(ECG_header) )
        ratio = 0;
        return
    end

    if( nargin < 3 || isempty(start_end) )
        start_end = [1 ECG_header.nsamp];
    end

    start_sample = start_end(1);
    end_sample = start_end(2);

    pb = progress_bar( 'Calculating ratios', 'Start' );

    % co_ocurrences respect other leads
    co_ocurrence = calc_co_ocurrences(time_serie);

    dummy = prdataset([]);
    dummy = prmapping([]);
    
    % load the trained classifier.
    aux_load = load('tc_qrs_q.mat');

    w_lablist = cellstr(getlabels(aux_load.wTrained_Classifier));
    FN_lab = find(strcmpi(w_lablist, 'FN'));
    FP_lab = find(strcmpi(w_lablist, 'FP'));
    TP_lab = find(strcmpi(w_lablist, 'TP'));

    lreferences = length(time_serie);

    pb.Loop2do = lreferences;
    pb.checkpoint('Calculating ratios.');
    
    ratio = nan(lreferences,1);

    estimated_labs = cell(lreferences,1);

    for ii = 1:lreferences

        pb.start_loop();

        this_ts = colvec(time_serie{ii});

        if( length(this_ts) < (ECG_header.nsamp / max_pattern_separation) )
%             warning( [mfilename ':few_anns'], 'Few detections in  %d\n', ii );
            ratio(ii) = 0;
        else
            
            pb.checkpoint([]);
            
            k_gaps = calc_k_gaps(this_ts);

            % build the feature matrix [ RR_i-1 RR_i RR_10 RR_60 co_ocurrences]
            % RR_i-1
            this_ts = round(this_ts * 1000 / ECG_header.freq);

            pb.checkpoint([]);
            
            aux_fm = [ calculate_RR_features(this_ts) colvec(co_ocurrence{ii}) ];

            pb.checkpoint([]);
            
            ds_result = prdataset(aux_fm) * aux_load.wTrained_Classifier;

            estimated_labs{ii} = renumlab(ds_result * labeld, char(w_lablist) );
            

            aux_val = sum(estimated_labs{ii} == TP_lab | estimated_labs{ii} == FN_lab);
            if( aux_val == 0 )
                this_se = 0;
            else
                this_se = sum(estimated_labs{ii} == TP_lab) / aux_val;
            end
            
            aux_val = sum(estimated_labs{ii} == TP_lab | estimated_labs{ii} == FP_lab);
            if( aux_val == 0 )
                this_pp = 0;
            else
                this_pp = sum(estimated_labs{ii} == TP_lab) / aux_val;
            end

            this_q = (2*this_se + this_pp)/3;
            ratio(ii) = k_gaps * this_q;
        end

        pb.end_loop();
    end

    clear pb
    
    function k = calc_k_gaps( time_serie )

        time_serie = time_serie( time_serie >= start_sample & time_serie <= end_sample );
        
        % add some samples in order to calculate the gaps among heartbeats (FN)
        if( time_serie(1) ~= start_sample )
            time_serie = [start_sample; time_serie];
        end

        if( time_serie(end) ~= end_sample )
            time_serie = [time_serie; end_sample];
        end

        time_serie = round(time_serie * 1000 / ECG_header.freq);
        
        RRserie = colvec(diff(time_serie));
        RRserie = [0;RRserie];

        gap_idx = find(RRserie > max_pattern_separation);
%         gap_end_idx = gap_idx;
%         gap_start_idx = gap_start_idx - 1;

        % Percent of the time we have gaps.
        k = 1 - (sum(RRserie(gap_idx)) / time_serie(end));

    end
    
    function fm = calculate_RR_features( time_serie )

        RRserie = diff(time_serie);
        RRserie = [RRserie(1); RRserie];
        
        % build the feature matrix [ RR_i-1 RR_i RR_10 RR_60 ]
        % RR_i-1
        fm = [ RRserie(1); RRserie(1:end-1) ];
        fm = [fm RRserie CalcFeatureRRx(time_serie, RRserie, 10000) CalcFeatureRRx(time_serie, RRserie, 60000)];

    end

    function start_end_aux = findStartEnd( bAux )

        start_aux = find(  bAux, 1, 'first' );
        end_aux = find(  bAux, 1, 'last' );
        start_end_aux = [start_aux end_aux];

    end

    function aux_mean = CalcFeatureRRx(time_serie, RRserie, win_size)
%         win_size in milliseconds

        aux_seq = 1:length(RRserie);
        aux_idx = arrayfun(@(a)( findStartEnd(  time_serie >= (time_serie(a) - win_size) & ... 
                                                time_serie <= time_serie(a) )), ...
                           aux_seq, 'UniformOutput', false);

        aux_mean = colvec(cellfun(@(a)(round(median(RRserie(a(1):a(2))))), aux_idx));
        
    end
    

end
