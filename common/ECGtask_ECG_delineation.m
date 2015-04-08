classdef ECGtask_ECG_delineation < ECGtask

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
        name = 'ECG_delineation';
        target_units = 'ADCu';
        doPayload = true;
        cAnnotationFields = { 'Pon' 'P' 'Poff' 'QRSon' 'qrs' 'Q' 'R' 'S' 'QRSoff' 'Ton' 'T' 'Toff' };
        
    end

    properties( GetAccess = public, SetAccess = private)
        % if user = memory;
        % memory_constant is the fraction respect to user.MaxPossibleArrayBytes
        % which determines the maximum input data size.
        memory_constant = 0.3;
        
        started = false;
        
    end
    
    properties( Access = private, Constant)
        
        cQRSdelineators = {'all-delineators' 'wavedet' };
        
    end
    
    properties( Access = private )
        
        delineators2do
        bWFDBdelineators
        tmp_path_local
        command_sep
        WFDB_bin_path
        lead_idx = [];
        lead_names
        
    end
    
    properties
        progress_handle
        only_ECG_leads = false;
        delineators = 'all-delineators';
        wavedet_config
        payload
        tmp_path
    end
    
    methods
           
        function obj = ECGtask_ECG_delineation(obj)

            % wavedet configuration structure
            obj.wavedet_config.setup.wavedet.QRS_detection_only = false;
            
        end
        
        function Start(obj, ECG_header, ECG_annotations)

            if( obj.only_ECG_leads )
                obj.lead_idx = get_ECG_idx_from_header(ECG_header);
            else
%                 'all-leads'
                obj.lead_idx = 1:ECG_header.nsig;
            end
            
            % lead names desambiguation
            str_aux = regexprep(cellstr(ECG_header.desc), '\W*(\w+)\W*', '$1');
            obj.lead_names = regexprep(str_aux, '\W', '_');

            [str_aux2, ~ , aux_idx] = unique( obj.lead_names );
            aux_val = length(str_aux2);
            
            if( aux_val ~= ECG_header.nsig )
                for ii = 1:aux_val
                    bAux = aux_idx==ii;
                    aux_matches = sum(bAux);
                    if( sum(bAux) > 1 )
                        obj.lead_names(bAux) = strcat( obj.lead_names(bAux), repmat({'v'}, aux_matches,1), cellstr(num2str((1:aux_matches)')) );
                    end
                end
            end
            
            obj.lead_names = regexprep(obj.lead_names, '\W*(\w+)\W*', '$1');
            obj.lead_names = regexprep(obj.lead_names, '\W', '_');
                        
            if( strcmpi('all-delineators', obj.delineators) )
                obj.delineators2do = obj.cQRSdelineators(2:end);
            else
                if( ischar(obj.delineators) )
                    obj.delineators2do = {obj.delineators};
                else
                    obj.delineators2do = obj.delineators;
                end
            end

            % local path required to avoid network bottlenecks in distributed filesystems 
            if( isunix() && exist('/scratch/', 'dir') )
                str_username = getenv('USER');
                obj.tmp_path_local = ['/scratch/' str_username filesep];
                if( ~exist(obj.tmp_path_local, 'dir') )
                    if(~mkdir(obj.tmp_path_local))
                        obj.tmp_path_local = '/scratch/';
                    end
                end
                obj.tmp_path = tempdir;
            else
                if( isempty(obj.tmp_path) )
                    obj.tmp_path = tempdir;
                    obj.tmp_path_local = tempdir;
                end
            end
            
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

                fnames = fieldnames(obj.payload);
                aux_idx = find(cell2mat( cellfun(@(a)(~isempty(strfind(a, 'corrected_'))), fnames, 'UniformOutput', false)));
                if( isempty(aux_idx) )

                    if( isfield(obj.payload, 'series_quality') )
                        [~, aux_idx] = sort(obj.payload.series_quality.ratios, 'descend');
                        aux_val = obj.payload.(obj.payload.series_quality.AnnNames{aux_idx(1),1}).(obj.payload.series_quality.AnnNames{aux_idx(1),2}) - ECG_start_offset + 1;
                        aux_val = aux_val( aux_val >= ECG_sample_start_end_idx(1) & aux_val < ECG_sample_start_end_idx(2) );
                        obj.payload = aux_val;
                        
                    else
                        for fname = rowvec(fnames)
                            if( isfield(obj.payload.(fname{1}), 'time') )
                                aux_val = obj.payload.(fname{1}).time - ECG_start_offset + 1;
                                aux_val = aux_val( aux_val >= ECG_sample_start_end_idx(1) & aux_val < ECG_sample_start_end_idx(2) );
                                obj.payload = aux_val;
                                break
                            end
                        end
                    end
                else
                    
                    for ii = rowvec(aux_idx)
                        aux_val = obj.payload.(fnames{ii}).time - ECG_start_offset + 1;
                        aux_val = aux_val( aux_val >= ECG_sample_start_end_idx(1) & aux_val < ECG_sample_start_end_idx(2) );
                        obj.payload = aux_val;
                    end
                end

            end
            
            cant_QRSdelineators = length(obj.delineators2do);
            
            for ii = 1:cant_QRSdelineators

                this_delineator = obj.delineators2do{ii};

                [this_delineator, this_delineator_name] = strtok(this_delineator, ':');
                if( isempty(this_delineator_name) )
                    this_delineator_name = this_delineator;
                else
                    this_delineator_name = this_delineator_name(2:end);
                end
                
                cprintf( 'Blue', [ 'Processing ECG delineator ' this_delineator_name '\n' ] );

                %% perform ECG delineation

                switch( this_delineator )

                    case 'wavedet'
                    %% Wavedet delineation

                        obj.progress_handle.checkpoint(['Processing ' this_delineator])
                    
                        try

                            [position_multilead, positions_single_lead] = wavedet_interface(ECG, ECG_header, obj.payload, obj.lead_idx, obj.wavedet_config, ECG_sample_start_end_idx, ECG_start_offset, obj.progress_handle);

                            for jj = 1:length(obj.lead_idx)
                                payload_out.wavedet.(obj.lead_names{obj.lead_idx(jj)}) = positions_single_lead(jj);
                            end
                            
                            payload_out.wavedet.multilead = position_multilead;

                        catch aux_ME

                            strAux = sprintf('Wavedet failed in recording %s\n', ECG_header.recname);
                            strAux = sprintf('%s\n', strAux);

                            report = getReport(aux_ME);
                            fprintf(2, '%s\nError report:\n%s', strAux, report);

                        end

                    case 'user'
                        %% user-defined delineator

                        obj.progress_handle.checkpoint(['Processing ' this_delineator_name])
                        
                        try

                            if( exist(this_delineator_name) == 2 )
                                
%                                 ud_func_pointer = eval(['@' this_delineator_name]);
                                ud_func_pointer = str2func(this_delineator_name);

                                obj.progress_handle.checkpoint([ 'User defined function: ' this_delineator_name])

                                ECG_header_aux = trim_ECG_header(ECG_header, obj.lead_idx);
                                
                                [positions_single_lead, position_multilead] = ud_func_pointer( double(ECG(:,obj.lead_idx)), ECG_header_aux, obj.progress_handle, obj.payload);

                                % filter and offset delineation
                                for jj = 1:length(obj.lead_idx)
                                    
                                    % filter heartbeats within range
                                    bAux = positions_single_lead(jj).qrs >= ECG_sample_start_end_idx(1) & positions_single_lead(jj).qrs <= ECG_sample_start_end_idx(2);

                                    aux_struct = [];
                                    for fn = rowvec(fieldnames(positions_single_lead(jj)))
                                        aux_val = positions_single_lead(jj).(fn{1});
                                        if( any(strcmpi(obj.cAnnotationFields, fn{1})) ) 
                                            aux_struct.(fn{1}) = aux_val(bAux) + ECG_start_offset - 1;
                                        else
                                            aux_struct.(fn{1}) = aux_val;
                                        end
                                    end
                                    
                                    payload_out.(this_delineator_name).(obj.lead_names{obj.lead_idx(jj)}) = aux_struct;
                                    
                                end
                                
                                if( ~isempty(position_multilead) )
                                    
                                    % filter heartbeats within range
                                    bAux = position_multilead.qrs >= ECG_sample_start_end_idx(1) & position_multilead.qrs <= ECG_sample_start_end_idx(2);

                                    aux_struct = [];
                                    for fn = rowvec(fieldnames(position_multilead))
                                        aux_val = position_multilead.(fn{1});
                                        if( any(strcmpi(obj.cAnnotationFields, fn{1})) ) 
                                            aux_struct.(fn{1}) = aux_val(bAux) + ECG_start_offset - 1;
                                        else
                                            aux_struct.(fn{1}) = aux_val;
                                        end
                                    end
                                    
                                    position_multilead = aux_struct;
                                    
                                    payload_out.(this_delineator_name).multilead = position_multilead;
                                    
                                end

                            else
                                disp_string_framed(2, sprintf('Function "%s" is not reachable in path.', this_delineator_name));
                                fprintf(1, 'Make sure that exist(%s) == 2\n',this_delineator_name);
                            end
                            
                        catch aux_ME

                            disp_string_framed(2, sprintf('Delineator "%s" failed in recording %s lead %s', this_delineator_name, ECG_header.recname, ECG_header.desc(jj,:) ) );                                

                            report = getReport(aux_ME);
                            fprintf(2, 'Error report:\n%s', report);

                        end

                end
    
            end
            
            
            
            
        end
        
        function payload = Finish(obj, payload, ECG_header)

            
        end
        
        function payload = Concatenate(obj, plA, plB)

            if( isempty(plA) )
                
                payload = plB;
                
            else
                for this_ECG_delineator = rowvec(fieldnames(plA))
                    
                    this_ECG_delineator = this_ECG_delineator{1};
                    
                    if( isfield(plB, this_ECG_delineator)  )
                        
                        for this_lead = rowvec(fieldnames(plA.(this_ECG_delineator)))
                            
                            this_lead = this_lead{1};
                        
                            if( isfield(plB.(this_ECG_delineator), this_lead)  )
                            
                                for fn = rowvec(fieldnames(plA.(this_ECG_delineator).(this_lead) ))

                                    if( isfield(plB.(this_ECG_delineator).(this_lead), fn{1}) )
                                        payload.( this_ECG_delineator).(this_lead).(fn{1}) = [ colvec(plA.( this_ECG_delineator).(this_lead).(fn{1}) ); colvec( plB.(this_ECG_delineator).(this_lead).(fn{1}) ) ];
                                    else
                                        error('ECGtask_ECG_delineation:BadTMPfiles', ['Results from ' this_ECG_delineator '.' this_lead  '.' fn{1} ' not found. Skipping concatenation.'])
                                    end

                                end
                            else
                                error('ECGtask_ECG_delineation:BadTMPfiles', ['Results from ' this_ECG_delineator '.' this_lead ' not found. Skipping concatenation.'] )
                            end
                        end
                        
                    else
                        error('ECGtask_ECG_delineation:BadTMPfiles', ['Results from ' this_ECG_delineator ' not found. Skipping concatenation.'] )
                    end
                end
            end
            

        end
        
        %% property restriction functions
        
        function set.delineators(obj,x)
            if( (ischar(x) || iscellstr(x)) )
                x = cellstr(x);
                aux_val = colvec(intersect(obj.cQRSdelineators, x));
                aux_idx = find(cellfun(@(a)(~isempty(strfind(a, 'user:'))), x));
                obj.delineators = [aux_val; colvec(x(aux_idx)) ];
            else
                warning('ECGtask_ECG_delineation:BadArg', 'Invalid delineators.');
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
                        warning('ECGtask_ECG_delineation:BadArg', ['Could not create ' x ]);
                    end
                end
                
            else
                warning('ECGtask_ECG_delineation:BadArg', 'tmp_path_local must be a string.');
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
                        warning('ECGtask_ECG_delineation:BadArg', ['Could not create ' x ]);
                    end
                end
                
            else
                warning('ECGtask_ECG_delineation:BadArg', 'tmp_path_local must be a string.');
            end
        end
        
    end
    
    methods ( Access = private )
        
        
    end
    
end
