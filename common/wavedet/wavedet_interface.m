function [multilead_positions single_lead_positions rhythm_parameters] = wavedet_interface(ECG_signal, ECG_header, beat_information, lead_to_delineate, wavedet_config, start_end, start_offset, progress_bar_hdl)

    if( nargin < 3 || ~(isstruct(beat_information) && isfield(beat_information, 'time')) )
        beat_information = [];
    else
        beat_information = beat_information.time;
    end

    if( nargin < 4 || isempty(lead_to_delineate))
        lead_to_delineate = 1:ECG_header.nsig;
    end

    if( nargin < 5 )
        wavedet_config = [];
    end
    
    if( nargin < 6 || isempty(start_end) )
        start_end = [1 size(ECG_signal,1) ];
    end
    
    if( nargin < 7 || isempty(start_offset) )
        start_offset = 1;
    end
    
    bProgress = true;
    if( nargin < 8 || isempty(progress_bar_hdl) )
        progress_bar_hdl = [];
        bProgress = false;
    end
    
    if( isfield(ECG_header, 'desc') )
        lead_names = ECG_header.desc;
    else
        lead_names = num2str((1:ECG_header.nsig)');
    end
    
    ecgkit_wave_defs();

    rhythm_parameters = [];
    multilead_positions = [];
    single_lead_positions = [];
    lead_processed_ok_idx = [];
    
    marks = cell(1,ECG_header.nsig);
%        corresponding respectively to P onset, P peak, P end, QRS onset, QRS
%        main peak, QRS end, T onset, T peak, T prima, T end

    % empty output struct
    for fn = cAnnotationOutputFields
        single_lead_positions.(fn{1}) = [];
    end
    single_lead_positions = repmat(single_lead_positions, ECG_header.nsig, 1);

    lead_processed_ok = 0;
    
    for ii = rowvec(lead_to_delineate)

%         fprintf(1, [ 'Processing lead ' lead_names(ii,:) '\n'])
%         fprintf(1, '.');
        
        try
            
            if(bProgress)
                progress_bar_hdl.checkpoint(['Lead ' lead_names(ii,:)])
            else
                fprintf(1, '.');
            end
            
            aux_struct = wavedet_3D( double(ECG_signal(:,ii)), beat_information, ECG_header, wavedet_config);
        
        catch ME
            
            report = getReport(ME);
            fprintf(2, '%s\n\n', report);
            fprintf(2, [ 'Problem found in lead ' lead_names(ii,:) '. Skipping this lead.\n'])
            
            continue
            
        end

%         % delete fields with no annotations
%         for field_nm = rowvec(fieldnames(aux_struct))
%             if( all(isnan(aux_struct.(field_nm{1}))) )
%                 aux_struct = rmfield(aux_struct, field_nm);
%             end
%         end
        
        % this field is necessary for the rules used later
        if( ~isfield( aux_struct, 'Tprima') && isfield( aux_struct, 'T') )
            aux_struct.Tprima = aux_struct.T;
        end

        aux_marks = positions2matrix(aux_struct, cAnnotationSLRfields) + start_offset - 1;
        
        % filter heartbeats within range
        bAux = aux_struct.qrs >= start_end(1) & aux_struct.qrs <= start_end(2);
        
        for fn = cAnnotationOutputFields
            if( isfield(aux_struct, fn{1}) ) 
                aux_val = aux_struct.(fn{1});
                aux_struct2.(fn{1}) = aux_val(bAux) + start_offset - 1;
            else
                aux_struct2.(fn{1}) = [];
            end
        end
        
        lead_processed_ok = lead_processed_ok + 1;
    
        single_lead_positions(ii) = aux_struct2;
        marks{ii} = aux_marks;
        lead_processed_ok_idx = [lead_processed_ok_idx; ii];
        
    end

    if( lead_processed_ok == 0 )
        warning('wavedet_cardiolund:DelineationImpossible', 'Delineation not possible.\n')
        return
    end
        
    if( lead_processed_ok >= 3 )
        
        if(bProgress)
            progress_bar_hdl.checkpoint(['Lead ' lead_names(ii,:)])
        else
            fprintf(1, '.\n');
        end
        
%         fprintf(1, 'Processing multilead rules\n');
        
        %Multilead rules.

        SLRmarks = SLR(marks(lead_processed_ok_idx),ECG_header.freq);

        multilead_positions.Pon = round(SLRmarks(:,1));
        multilead_positions.P = round(SLRmarks(:,2));
        multilead_positions.Poff = round(SLRmarks(:,3));
        multilead_positions.QRSon = round(SLRmarks(:,4));
        multilead_positions.qrs = round(SLRmarks(:,5));
        multilead_positions.QRSoff = round(SLRmarks(:,6));
        multilead_positions.Ton = round(SLRmarks(:,7));
        multilead_positions.T = round(SLRmarks(:,8));
        multilead_positions.Tprima = round(SLRmarks(:,9));
        multilead_positions.Toff = round(SLRmarks(:,10));

        % trust more on peaks.
        multilead_positions.Ton( multilead_positions.Ton >= multilead_positions.T ) = nan;
        multilead_positions.Toff( multilead_positions.Toff <= multilead_positions.T ) = nan;
        multilead_positions.Pon( multilead_positions.Pon >= multilead_positions.P ) = nan;
        multilead_positions.Poff( multilead_positions.Poff <= multilead_positions.P ) = nan;
        multilead_positions.QRSon( multilead_positions.QRSon >= multilead_positions.qrs ) = nan;
        multilead_positions.QRSoff( multilead_positions.QRSoff <= multilead_positions.qrs ) = nan;
        
    else
        
        % in case user want multilead, and not possible.
        if( length(lead_to_delineate) >= 3 )
            fprintf(2, [ 'Delineation not possible in several leads. Using the delineation of lead ' strtrim(lead_names(lead_processed_ok_idx(1),:)) '. Consider reviewing the delineation.\n'])
        else
            fprintf(1, [ 'No multilead strategy used, instead using the delineation of lead ' strtrim(lead_names(lead_processed_ok_idx(1),:)) '.\n'])
        end
            
        multilead_positions = single_lead_positions(lead_processed_ok_idx(1));
        
    end
    
    if( nargout > 2 )
        % calculate some intervals
        if(length(multilead_positions.qrs) > 1)
            rhythm_parameters.RR = diff( multilead_positions.qrs );
            rhythm_parameters.RR = [rhythm_parameters.RR(1); rhythm_parameters.RR] / ECG_header.freq;
        else
            rhythm_parameters.RR = [];
        end

        rhythm_parameters.QT = (multilead_positions.Toff - multilead_positions.QRSon)/ ECG_header.freq;
        rhythm_parameters.QTp = (multilead_positions.T - multilead_positions.QRSon)/ ECG_header.freq;
    end
    