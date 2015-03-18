classdef ECGtask_PPG_ABP_corrector < ECGtask

% ECGtask for ECGwrapper (for Matlab)
% ---------------------------------
% 
% Description:
% 
% Abstract class for defining ECGtask interface
% 
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Birthdate  : 18/2/2013
% Last update: 18/2/2013
       
    properties(GetAccess = public, Constant)
        name = 'PPG_ABP_corrector';
        % Require the parent Wrapper object to this task
        target_units = 'Wrapper';
        doPayload = true;
    end

    properties( GetAccess = public, SetAccess = private)
        % if user = memory;
        % memory_constant is the fraction respect to user.MaxPossibleArrayBytes
        % which determines the maximum input data size.
        % Size > 1 means that this task can handle big data structures, as
        % this case.       
        memory_constant = realmax;
        
        started = false;
        
    end
    
    properties( Access = private, Constant)
        
        fig_hdl = 1;
        cAnnotationsFieldNamesRequired = {'time', 'qrs' };
        
    end
    
    properties( Access = private )
        
    end
    
    properties
        
        progress_handle
        caller_variable = 'payload'
        tmp_path
        payload
        signal_payload = false
        
    end
    
    methods
           
        function obj = ECGtask_PPG_ABP_corrector(obj)
            
        end
        
        function Start(obj, ECG_header, ECG_annotations)
            
            obj.started = true;

        end
        
        function payload = Process(obj, ECG, ECG_start_offset, ECG_sample_start_end_idx, ECG_header, ECG_annotations, ECG_annotations_start_end_idx )

            payload = [];

            obj.progress_handle.hide()
            
%             if( ~obj.started )
%                 obj.Start(ECG_header);
%                 if( ~obj.started )
%                     cprintf('*[1,0.5,0]', 'Task %s unable to be started for %s.\n', obj.name, ECG_header.recname);
%                     return
%                 end
%             end
            
            this_start_end = (ECG_sample_start_end_idx + ECG_start_offset - 1);
%             aux_val = this_start_end ./ ECG_header.freq;
%             disp_string_title(1, sprintf( 'Correcting from %s to %s', Seconds2HMS(aux_val(1)), Seconds2HMS(aux_val(2)) ) );
            
            if( isempty(obj.payload) )
                Ann_struct = ECG_annotations;
            else
                Ann_struct = [];
                
                if( ~isempty(obj.payload) )
                    
                    if( isfield(obj.payload, 'series_quality' ) )
                        Ann_struct.series_quality = obj.payload.series_quality;
                    end
                    
                    % filter and move QRS detections in the input payload
                    for fn = rowvec(fieldnames(obj.payload))
                        if( isstruct( obj.payload.(fn{1}) ) )
                            aux_val = intersect( fieldnames(obj.payload.(fn{1})), obj.cAnnotationsFieldNamesRequired );
                            if( ~isempty(aux_val) )
                                aux_val = obj.payload.(fn{1}).(aux_val{1});
                                aux_val2 = size(aux_val);

                                if( aux_val2(1) == 1 && aux_val2(2) > 1 )
                                    aux_val = colvec(aux_val);
                                elseif( aux_val2(1) > 1 && aux_val2(2) > 1 && aux_val2(2) > aux_val2(1) )
                                    aux_val = aux_val';
                                end

                                if( size(aux_val, 2) == 1 )
                                    aux_val = aux_val(aux_val >= this_start_end(1) & aux_val <= this_start_end(2)) - ECG_start_offset + 1 ;

                                elseif( size(aux_val, 2) > 1 )
                                    bAux = aux_val(:,1) >= this_start_end(1) & aux_val(:,1) <= this_start_end(2);
                                    aux_val = aux_val(bAux, :);
                                    aux_val(:,1) = aux_val(:, 1) - ECG_start_offset + 1 ;
                                end
                                Ann_struct.(fn{1}).time = aux_val;
                            end
                        end
                    end
                end
                
            end
            
            if( ~isempty(ECG_start_offset) || ECG_start_offset > 1)
                 [~, iHours, iMins, iSeconds, iMilli ] = Seconds2HMS(ECG_start_offset/ECG_header.freq);
                 ECG_header.btime = datestr([2014,1,1,iHours, iMins,iSeconds], 'HH:MM:SS');
            end

            QRScorrector('ECG', ECG, 'QRS_annotations', Ann_struct, 'Figure', figure(obj.fig_hdl) );

            disp_string_framed('*Blue', 'User interaction required' );

            aux_str = ['<a href="matlab:figure(' num2str(obj.fig_hdl) ')">figure ' num2str(obj.fig_hdl) '</a>'];
            
            aux_str2 = '<a href="matlab:dbcont">F5 (Run)</a>';
            
            fprintf(1, 'This ECGtask allow user interaction. Press [CTRL + G] in %s to save results and press %s to continue.\n', aux_str, aux_str2)
            keyboard

            % last chance to save results
            if( ishandle(obj.fig_hdl) && (isempty(payload) || ~isstruct(payload)) )
                disp_string_framed('*[1,0.5,0]', 'Payload variable not saved' );
                fprintf(1, 'Press [CTRL + G] in %s to save results and press %s to continue.\n', aux_str, aux_str2)
                keyboard

                if( ~isempty(payload) && isstruct(payload) )
                    % save work done, to further improve in following
                    % invocations.
                    disp_string_framed('*Magenta', 'Corrections saved' );

                else
                    disp_string_framed('*Red', 'Data not saved' );
                end   
                
            else
                disp_string_framed('*Magenta', 'Corrections saved' );
            end

            if(ishandle(obj.fig_hdl))
                close(obj.fig_hdl)
            end
            
            if(~(isempty(payload) || ~isstruct(payload)))
                % move QRS corrections according to the input offset
                for fn = rowvec(fieldnames(payload))
                    if( isfield(payload.(fn{1}), 'time' ) )
                        payload.(fn{1}).time = payload.(fn{1}).time + ECG_start_offset - 1;
                    end
                end
            end
            
            obj.progress_handle.show()
            
        end
        
        function payload = Finish(obj, payload, ECG_header)
           
        end
        
        function payload = Concatenate(obj, plA, plB)

            if( isempty(plA) )
                
                payload = plB;
                
            else
                
                fields = rowvec(unique( [rowvec(fieldnames(plA)) rowvec(fieldnames(plB))] ));
                
                diff_fields = setdiff(fields, {'series_quality'});
                
                for fn = diff_fields
                    
                    if( isfield(plA, fn{1}) && isfield(plB, fn{1}) )
                        payload.(fn{1}).time = [ colvec(plA.(fn{1}).time); colvec(plB.(fn{1}).time) ];
                    
                    elseif( ~isfield(plA, fn{1}) && isfield(plB, fn{1}) )
                        payload.(fn{1}).time = colvec(plB.(fn{1}).time);
                    end
                    
                end
                
                aux_idx = find(strcmpi(fields, 'series_quality'));
                
                if( ~isempty(aux_idx) )
                    
                    if( isfield(plA, 'series_quality') && isfield(plB, 'series_quality') )
                        payload.series_quality.ratios = [plA.series_quality.ratios plB.series_quality.ratios];
                        payload.series_quality.estimated_labs = cellfun(@(a,b)( [colvec(a);colvec(b)] ) , plA.series_quality.estimated_labs, plB.series_quality.estimated_labs, 'UniformOutput', false);
                        payload.series_quality.AnnNames = plA.series_quality.AnnNames;
                    end

                end
                
            end

        end

        %% property restriction functions

        
    end
    
    methods ( Access = private )
        
        
    end
    
end
