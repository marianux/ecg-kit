function [multilead_positions single_lead_positions rhythm_parameters] = wavedet_interface(ECG_signal, ECG_header, beat_information, lead_to_delineate, wavedet_config)

    if( nargin < 3 )
        beat_information = [];
    end

    if( nargin < 4 || isempty(lead_to_delineate))
        lead_to_delineate = 1:ECG_header.nsig;
    end

    if( nargin < 5 )
        wavedet_config = [];
    end
    
    if( isfield(ECG_header, 'desc') )
        lead_names = ECG_header.desc;
    else
        lead_names = num2str((1:ECG_header.nsig)');
    end
    
    cAnnotationFieldNames = { 'Pon' 'P' 'Poff' 'QRSon' 'qrs' 'QRSoff' 'Ton' 'T' 'Tprima' 'Toff' };

    rhythm_parameters = [];
    multilead_positions = [];
    single_lead_positions = [];
    
    marks = cell(1,ECG_header.nsig);
%        corresponding respectively to P onset, P peak, P end, QRS onset, QRS
%        main peak, QRS end, T onset, T peak, T prima, T end

    lead_processed_ok = 0;
    
    for ii = lead_to_delineate

%         fprintf(1, [ 'Processing lead ' lead_names(ii,:) '\n'])
        fprintf(1, '.');
        
        try
            
            aux_struct = wavedet_3D(ECG_signal(:,ii), beat_information, ECG_header, wavedet_config);
        
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

        % get the field with max num of annotations
        aux_maxlength = 0;
        for fname = rowvec(fieldnames(aux_struct))
            aux_maxlength = max(aux_maxlength, length(aux_struct.(fname{1})));
        end
        
        %store to use rules on all annotations at the end
        aux_marks = nan(aux_maxlength,10);
        
        count = 1;
        for fn = cAnnotationFieldNames
            fn = fn{1};
            if( isfield( aux_struct, fn) )
                aux_marks(1:length(aux_struct.(fn)),count) = colvec(aux_struct.(fn));
            end
            count = count + 1;
        end

        for fn = cAnnotationFieldNames
            if( isfield(aux_struct, fn{1}) ) 
                aux_struct2.(fn{1}) = aux_struct.(fn{1});
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
        
%         fprintf(1, 'Processing multilead rules\n');
        fprintf(1, '.\n');
        
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
            fprintf(2, [ 'To few leads to use multilead rules. Using the delineation of lead ' lead_names(1,:) '. Consider reviewing the delineation.\n'])
        else
            fprintf(1, [ 'Using the delineation of lead ' lead_names(lead_processed_ok_idx(1),:) '.\n'])
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
    