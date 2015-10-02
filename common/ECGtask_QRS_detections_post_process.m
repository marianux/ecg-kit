classdef ECGtask_QRS_detections_post_process < ECGtask

% ECGtask for ECGwrapper (for Matlab)
% ---------------------------------
% 
% Description:
% 
% Abstract class for defining ECGtask interface
% 
% Adding user-defined QRS detectors:
% A QRS detector that has the following interface can be added to the task:
% 
%     [positions_single_lead, position_multilead] = your_ECG_delineator( ECG_matrix, ECG_header, progress_handle, payload_in);
% 
% where the arguments are:
%    + ECG_matrix, is a matrix size [ECG_header.nsamp ECG_header.nsig]
%    + ECG_header, is a struct with info about the ECG signal, such as:
%         .freq, the sampling frequency
%         .desc, description about the signals.
%    + progress_handle, is a handle to a waitbar object, that can be used
%          to track the progress within your function. See the
%          documentation about this class in this kit.
%    + payload_in, is a user data variable allowed to be sent each call to
%          your function. It is sent, via the payload property of this
%          class, for example: 
% 
%         this_ECG_wrappers.ECGtaskHandle.payload = your_variable;
%         this_ECG_wrappers.ECGtaskHandle.payload = {your_var1 your_var2};
%         this_ECG_wrappers.ECGtaskHandle.payload = load(cached_filenames);
% 
%          In the context of delineation, it is thought to be a user
%          corrected, or "gold quality" QRS location, in order to 
%          improve the wave delineation quality. If "payload_in" is a
%          struct, this function will automatically filter and time-shift
%          all QRS detection fields started with the string "corrected_".
%          For this purpose, QRScorrector task, automatically appends this
%          string to eachm anually reviewed QRS location series. 
% 
% the output of your function must be:
%    + positions_single_lead, a cell array size ECG_header.nsig with the
%          QRS sample locations found in each lead.
%    + position_multilead, a numeric vector with the QRS locations
%          calculated using multilead rules.
% 
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Birthdate  : 18/2/2013
% Last update: 18/2/2013
       
    properties(GetAccess = public, Constant)
        name = 'QRS_detections_post_process';
        target_units = 'ADCu';
        doPayload = true;
    end

    properties( GetAccess = public, SetAccess = private)
        % if user = memory;
        % memory_constant is the fraction respect to user.MaxPossibleArrayBytes
        % which determines the maximum input data size.
        memory_constant = 0.3;
        
        started = false;
        
    end
    
    properties( Access = private, Constant)
        
    end
    
    properties( Access = private )
        
        tmp_path_local
                
    end
    
    properties
        post_proc_func
        progress_handle
        payload
        tmp_path
        CalculatePerformance = false;
        
    end
    
    methods
           
        function obj = ECGtask_QRS_detections_post_process (obj)
            
        end
        
        function Start(obj, ECG_header, ECG_annotations)
            
            obj.started = true;
            
        end
        
        function payload_out = Process(obj, ECG, ECG_start_offset, ECG_sample_start_end_idx, ECG_header, ECG_annotations, ECG_annotations_start_end_idx )
            
            payload_out = [];
            
            if( ~obj.started )
                obj.Start(ECG_header);
                if( ~obj.started )
                    cprintf('*[1,0.5,0]', 'Task %s unable to be started for %s.\n', obj.name, ECG_header.recname);
                    return
                end
            end
            
            % payload property is used in this task to input an external QRS
            % detector, or manually corrected detections.
            if( isstruct(obj.payload) )

                AnnNames = obj.payload.series_quality.AnnNames(:,1);
                aux_struct = obj.payload;
                
                for fn = rowvec(AnnNames)
                    aux_val = obj.payload.(fn{1}).time - ECG_start_offset + 1;
                    aux_val = aux_val( aux_val >= ECG_sample_start_end_idx(1) & aux_val < ECG_sample_start_end_idx(2) );
                    aux_struct.(fn{1}).time = aux_val;
                end

            else
                return
            end

            for this_func = rowvec(obj.post_proc_func)
                
                try

                    obj.progress_handle.checkpoint([ 'User defined function: ' this_func{1}])

    %                 this_func_ptr = eval(['@' this_func]);
                    this_func_ptr = str2func(this_func{1});

                    this_payload = this_func_ptr( aux_struct, ECG_header, ECG_sample_start_end_idx );

                    for fn = rowvec(fieldnames(this_payload))
                        payload_out.(fn{1}) = this_payload.(fn{1});
                    end
                    
                catch aux_ME

                    disp_string_framed(2, sprintf('User-function "%s" failed in recording %s', this_func{1}, ECG_header.recname ) );                                

                    report = getReport(aux_ME);
                    fprintf(2, 'Error report:\n%s', report);

                end
                
            end

            obj.progress_handle.checkpoint('Adding quality metrics')

            % Add QRS detections quality metrics, Names, etc.
            payload_out = calculateSeriesQuality(payload_out, ECG_header, ECG_sample_start_end_idx );

            % offset in time
            AnnNames = payload_out.series_quality.AnnNames(:,1);
            for fn = rowvec(AnnNames)
                aux_val = payload_out.(fn{1}).time + ECG_start_offset - 1;
                payload_out.(fn{1}).time = aux_val;
            end

            % calculate performance
            if( obj.CalculatePerformance )
                AnnNames = payload_out.series_quality.AnnNames(:,1);
                cant_lead_name = size(AnnNames,1);
                payload_out.series_performance.conf_mat = zeros(2,2,cant_lead_name);
                payload_out.series_performance.error = nan(cant_lead_name,2);

                if(isempty(ECG_annotations)) 
                    disp_string_framed(2, sprintf('Trusted references not found for %s', ECG_header.recname) );
                else
                    % offset refs, produced anns were already shifted
                    ECG_annotations.time = ECG_annotations.time + ECG_start_offset - 1;
                    for kk = 1:cant_lead_name
                        this_lead_det = payload_out.(AnnNames{kk}).time;
                        [payload_out.series_performance.conf_mat(:,:,kk), TP_gs_idx, TP_det_idx] = bxb(ECG_annotations, this_lead_det, ECG_header );
                        % meadian/mad of the error wrt the gold standard
                        aux_val = this_lead_det(TP_det_idx) - ECG_annotations(TP_gs_idx);
                        payload_out.series_performance.error(kk,:) = [ median(aux_val) mad(aux_val) ];
                    end            
                end

            end


            obj.progress_handle.checkpoint('Done')
            
        end
        
        function payload = Finish(obj, payload, ECG_header)

            if( isfield(payload, 'series_quality') && isfield(payload.series_quality, 'ratios') )
                payload.series_quality.ratios = mean(payload.series_quality.ratios, 2);
            end
            
        end
        
        function payload = Concatenate(obj, plA, plB)

            payload = ConcatenateQRSdetectionPayloads(obj, plA, plB);

        end
        
        %% property restriction functions

        function set.post_proc_func(obj, x)
            
            if( ischar(x) )
                
                if( exist(x) == 2 )
                    obj.post_proc_func = cellstr(x);
                else
                    disp_string_framed(2, sprintf('Function "%s" is not reachable in path.', x));  
                    fprintf(1, 'Make sure that exist(%s) == 2\n',x);
                end
                
            elseif( iscellstr(x) )

                if( any( cellfun( @(a)(exist(a)), x) == 2) )
                    obj.post_proc_func = x;
                else
                    disp_string_framed(2, sprintf('Function "%s" is not reachable in path.', x));  
                    fprintf(1, 'Make sure that exist(%s) == 2\n',x);
                end
                
            else
                warning('ECGtask_QRS_detections_post_process:BadArg', 'post_proc_func must be a string.');
            end
        end
        
                
        function set.tmp_path_local(obj,x)
            if( ischar(x) )
                if(exist(x, 'dir'))
                    obj.tmp_path_local = x;
                else
                    if(mkdir(x))
                        obj.tmp_path_local = x;
                    else
                        warning('ECGtask_QRS_detections_post_process:BadArg', ['Could not create ' x ]);
                    end
                end
                
            else
                warning('ECGtask_QRS_detections_post_process:BadArg', 'tmp_path_local must be a string.');
            end
        end
        
        function set.tmp_path(obj,x)
            if( ischar(x) )
                if(exist(x, 'dir'))
                    obj.tmp_path = x;
                else
                    if(mkdir(x))
                        obj.tmp_path = x;
                    else
                        warning('ECGtask_QRS_detections_post_process:BadArg', ['Could not create ' x ]);
                    end
                end
                
            else
                warning('ECGtask_QRS_detections_post_process:BadArg', 'tmp_path_local must be a string.');
            end
        end
        
    end
    
    methods ( Access = private )
        
        
    end
    
end
