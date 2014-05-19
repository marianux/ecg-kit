classdef ECGtask_PCA_proj_basis < ECGtask

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
        name = 'PCA_proj_basis';
        target_units = 'uV';
        doPayload = false;
    end

    properties( GetAccess = public, SetAccess = private)
        % if user = memory;
        % memory_constant is the fraction respect to user.MaxPossibleArrayBytes
        % which determines the maximum input data size.
        memory_constant = 0.3;
    end
    
    properties( Access = private, Constant)
    
        %sample 10 minutes from the whole recording.
        time_sample = 10 * 60; %seconds
        time_heartbeat_window = 2; %seconds around the heartbeat
    end
    
    properties( Access = private )
        ECG_slices
        QRS_sample_idx
        wavelet_filters
        some_window
        halfwin_samples
    end
    
    properties
        progress_handle 
        cant_QRS_locations
        autovec
        user_string = 'results'
    end
    
    methods
           
        function obj = ECGtask_PCA_proj_basis(obj)
            obj.ECG_slices = [];
        end
        
        function Start(obj, ECG_header, ECG_annotations)

            % ECG cell for grouping QRS complex slices
            obj.ECG_slices = [];
                                    
            % Sample which QRS complexes will consider into de PCA basis
            % calculation
            cant_QRS_sample = round(obj.time_sample / obj.time_heartbeat_window);

            obj.QRS_sample_idx = sort(randsample(obj.cant_QRS_locations, cant_QRS_sample));
            
            obj.halfwin_samples = round( obj.time_heartbeat_window/2*ECG_header.freq);
            
            obj.some_window = colvec(blackman(2*obj.halfwin_samples+1));
            
        end
        
        function payload = Process(obj, ECG, ECG_sample_start_end_idx, ECG_header, ECG_annotations, ECG_annotations_start_end_idx )
            
            % this object doesn't generate any payload
            payload = [];

            this_iter_QRS_sample_idx = find(obj.QRS_sample_idx >= ECG_annotations_start_end_idx(1) & obj.QRS_sample_idx <= ECG_annotations_start_end_idx(2));
            
            if( ~isempty(this_iter_QRS_sample_idx) )
                
                % change the reference to this iteration
                this_iter_QRS_sample_idx = obj.QRS_sample_idx(this_iter_QRS_sample_idx) - ECG_annotations_start_end_idx(1) + 1;
                
                QRS_locations = ECG_annotations.time;

                this_iter_ECG_size = size(ECG,1);

                ECG = BaselineWanderRemovalSplines( ECG, QRS_locations, ECG_header.freq);

                aux_idx = arrayfun(@(a)( max(1, QRS_locations(a) - obj.halfwin_samples): ...
                                         min( this_iter_ECG_size, QRS_locations(a) + obj.halfwin_samples)) , ...
                                   colvec(this_iter_QRS_sample_idx), 'UniformOutput', false);

                obj.ECG_slices = [obj.ECG_slices; cellfun(@(a)( bsxfun(@times, ECG(a,:), obj.some_window) ), aux_idx, 'UniformOutput', false)];
            
            end
            
        end
        
        function Finish(obj, ECG_header)
            
            ECG = cell2mat(obj.ECG_slices);

            % Wavelet transform calculation
            wtECG = qs_wt(ECG, 4, ECG_header.freq, obj.wavelet_filters);

            % calculate PCA in wt scale 4 of the ECG
            result = mcdcov( wtECG, 'plots', 0);
            WTecg_cov = result.cov;
            [obj.autovec autoval] = eig(WTecg_cov); 
            autoval = diag(autoval);
            [~, autoval_idx] = sort(autoval, 'descend');
            obj.autovec = obj.autovec(:,autoval_idx);   
            
        end
        
        function payload = Concatenate(obj, plA, plB)
            
            payload = [];
            % not implemented
            
        end

    end
    
end
