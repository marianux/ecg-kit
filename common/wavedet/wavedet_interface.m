function [multilead_positions single_lead_positions rhythm_parameters] = wavedet_interface(ECG_signal, ECG_header, beat_information, lead_to_delineate, wavedet_config, start_end, start_offset, progress_bar_hdl)

    if( nargin < 3 )
        beat_information = [];
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
    
    cAnnotationOutputFields = { 'Pon' 'P' 'Poff' 'Ptipo' 'QRSon' 'qrs' 'Q' 'R' 'S' 'QRSoff' 'Ton' 'T' 'Tprima' 'Toff' 'Ttipo' };
    cAnnotationSLRfields = { 'Pon' 'P' 'Poff' 'QRSon' 'qrs' 'QRSoff' 'Ton' 'T' 'Tprima' 'Toff' };

    rhythm_parameters = [];
    multilead_positions = [];
    single_lead_positions = [];
    
    marks = cell(1,ECG_header.nsig);
%        corresponding respectively to P onset, P peak, P end, QRS onset, QRS
%        main peak, QRS end, T onset, T peak, T prima, T end

    lead_processed_ok = 0;
    
    for ii = lead_to_delineate

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

        aux_marks = positions2matrix(aux_struct, cAnnotationSLRfields);
        
        % filter heartbeats within range
        bAux = aux_struct.qrs >= start_end(1) & aux_struct.qrs <= start_end(2);
        
        for fn = cAnnotationOutputFields
            if( isfield(aux_struct, fn{1}) ) 
                aux_val = aux_struct.(fn{1});
                aux_struct2.(fn{1}) = aux_val(bAux) - start_offset + 1;
            else
                aux_struct2.(fn{1}) = [];
            end
        end
        
        if( lead_processed_ok == 0 )
            single_lead_positions = aux_struct2;
            lead_processed_ok_idx = ii;
        else
            single_lead_positions(lead_processed_ok+1) = aux_struct2;
            lead_processed_ok_idx = [lead_processed_ok_idx; ii];
        end
        
        marks{lead_processed_ok+1} = aux_marks;
        
        lead_processed_ok = lead_processed_ok + 1;
        
    end

    if( lead_processed_ok == 0 )
        warning('wavedet_cardiolund:DelineationImpossible', [ 'Delineation not possible.\n'])
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

        SLRmarks = SLR(marks,ECG_header.freq);

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
        
    else
        
        % in case user want multilead, and not possible.
        if( length(lead_to_delineate) >= 3 )
            fprintf(2, [ 'Delineation not possible in several leads. Using the delineation of lead ' lead_names(lead_processed_ok_idx(1),:) '. Consider reviewing the delineation.\n'])
        else
            fprintf(1, [ 'No multilead strategy used, instead using the delineation of lead ' lead_names(lead_processed_ok_idx(1),:) '.\n'])
        end

        aux_marks = marks{1};

        multilead_positions.Pon = round(aux_marks(:,1));
        multilead_positions.P = round(aux_marks(:,2));
        multilead_positions.Poff = round(aux_marks(:,3));
        multilead_positions.QRSon = round(aux_marks(:,4));
        multilead_positions.qrs = round(aux_marks(:,5));
        multilead_positions.QRSoff = round(aux_marks(:,6));
        multilead_positions.Ton = round(aux_marks(:,7));
        multilead_positions.T = round(aux_marks(:,8));
        multilead_positions.Tprima = round(aux_marks(:,9));
        multilead_positions.Toff = round(aux_marks(:,10));        
            
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
    