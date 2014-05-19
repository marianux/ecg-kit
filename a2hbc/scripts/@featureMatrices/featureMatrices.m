classdef featureMatrices < handle

% featureMatrices object for ECGwrapper
% -------------------------------------
% 
% Description:
% 
% 
% Arguments: (specified as ECGwrapper('arg_name1', arg_val1, ... , 'arg_nameN', arg_valN) )
% 
%     a. ECG specification:
% 
%         a.1 Recording filename where ECG is stored in a valid format.
% 
%           + recording_name: ECG recording to be classified.
%           + recording_format : Valid ECG format. (MIT, ISHNE, AHA, HES, MAT)
% 
%     b. Operating modes
% 
%       
%     c. Modifiers
% 
% 
% Output:
% 
% Examples:
% 
% 
% Limits and Known bugs:
%   Probably a lot :( ... but dont panic! send me feedback if you need
%   help.
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Version: 0.1 beta
% Birthdate  : 20/3/2012
% Last update: 20/3/2012
       
    properties(GetAccess = private, Constant)
        sampling_rate = 360; %Hz
        fnPayload = {'featMat_clust' 'featMat_ldc' 'ECG' 'QRS_locations' 'RR_intervals'};
        fnFMfile = {'featMat_clust' 'featMat_ldc'};
        fnECGfile = {'ECG' 'QRS_locations' 'RR_intervals'};
        SmallValue = 0.0001;
    end
    
    properties ( Access = private )
        delayLP
        autovec
        q_filters
        LP
        LPsize
        scale_idx
        scales
        global_struct
    end
    
    properties (SetAccess = private, GetAccess = public)  
        %read-only 
        name
        class_labeling
        true_labels
        %auxiliary indexes
        file_count = 1;
        File_list = [];
        File_idx = [];
        QRS_idx = [];
        
    end

    properties
        %read-write
        recording_name
        recording_format
        progress_hdl
        tmp_path
        bCache
        cant_pids
        this_pid
        cant_iteration
        this_iteration
    end
    

    methods 
        
        function obj = featureMatrices(class_labeling)

            %% Constants and definitions

            obj.name = 'feat_matrix_construction';
            obj.class_labeling = class_labeling;

        end
        
        function Start(obj)
        

            %% Calculate accesories signals
            
            %load low-pass preprocessing filter.
            load( ['LP Fs ' num2str(obj.sampling_rate) 'Hz - Fc 35Hz - 80db stop.mat' ]);
            obj.LPsize = length(LP.numerator);
            obj.delayLP = round((obj.LPsize-1)/2);
            obj.LP = LP;
            clear LP
            
            % Wavelet transform filter bank design
            obj.scales = 1:6;
            CantScales = length(obj.scales);
            MaxScales = max(obj.scales);
            obj.scale_idx= nan(MaxScales,1);

            filters_cache_filename = ['wt_filters_' num2str(MaxScales) ' scales_' num2str(obj.sampling_rate) ' Hz.mat' ];
            if( exist(filters_cache_filename, 'file') )
                aux = load( filters_cache_filename );
                obj.q_filters = aux.q_filters;
                clear aux
            else
                obj.q_filters = qs_filter_design(MaxScales, obj.sampling_rate);
            end
            
            for ii = 1:CantScales
                %indice para saber que escala corresponde a cada columna de MyWavelets
                obj.scale_idx(obj.scales(ii)) = ii;
            end

            % Update point
            obj.progress_hdl.checkpoint('PCA projection');

%             obj.autovec = PCA_proj_basis(obj.recording_name, obj.recording_format, obj.q_filters);
            
            obj.autovec = [-0.696609733483207 -0.403719836632163 -0.593081084444745;-0.487887339191286 -0.339521516938205 0.804171053814316;-0.526023595928265 0.849550135686930 0.0395441965529324;];
            
            % data structure to keep data through all iterations.
            obj.global_struct = [];
            
        end
        
        function payload = Process(obj, ECG, this_header, this_iter_ECG_relative_start_end_idx, this_iter_annotations)

            obj.progress_hdl.checkpoint('Label parsing');
            
            [QRS_locations, obj.true_labels ] = Annotation_process(this_iter_annotations, obj.recording_format, obj.class_labeling);
            
            % RR interval sequence
            RR_intervals = diff(QRS_locations, 1);
            RR_intervals = colvec([ RR_intervals(1); RR_intervals ]);

            dRR_intervals = diff(RR_intervals, 1);
            dRR_intervals = colvec([ dRR_intervals(2); dRR_intervals(2); dRR_intervals(2:end) ]);

            % interval resampling.
            RR_intervals = RR_intervals * obj.sampling_rate / this_header.freq;
            dRR_intervals = dRR_intervals * obj.sampling_rate / this_header.freq;
            
            %% Preprocessing
            
            % Update point
            obj.progress_hdl.checkpoint('Resampling');
            
            cant_QRS_locations = length(QRS_locations);
            this_iter_QRS_seq_idx = 1:cant_QRS_locations;

            %resample to obj.sampling_rate Hz
            ECG = resample(ECG, obj.sampling_rate, this_header.freq);
            QRS_locations = colvec(round(QRS_locations * obj.sampling_rate / this_header.freq));
        %         this_iter_true_labels = true_labels(this_iter_QRS_seq_idx);
            this_iter_ECG_resampled_size = size(ECG,1);

            if( this_header.nsig > 2 ) 
                %multilead approach.
                % project ECG and wtECG to obtain PCA components
                ECG = ECG * obj.autovec(:,1:2);
        %             wtECG = cellfun(@(a)( a * obj.autovec(:,1:2) ), mat2cell(wtECG, this_iter_ECG_resampled_size, obj.ECG_header.nsig, ones(1,CantScales) ), 'UniformOutput', false);
        %             wtECG = cell2mat(wtECG);
            end

            % Update point
            obj.progress_hdl.checkpoint('Filtering');

            %Low pass filtering @ 35Hz => ECG recovery
            ECG = filter(obj.LP, ECG);
            %delay compensation.
            ECG = [ zeros(obj.delayLP, size(ECG,2) ) ;ECG(obj.LPsize:end,:) ; zeros(obj.delayLP, size(ECG,2))];

            %Quito la linea de base.
        %     ECG = BaselineWanderRemovalMedian( ECG, obj.sampling_rate);
            ECG = BaselineWanderRemovalSplines( ECG, QRS_locations, obj.sampling_rate);

            % Update point
            obj.progress_hdl.checkpoint('Wavelet transform calculation');

            wtECG = qs_wt(ECG, obj.scales, obj.sampling_rate, obj.q_filters);

            %% Features Calculation

            % Update point
            obj.progress_hdl.checkpoint('Features calculation');

            this_iter_seq_idx = 1:length(QRS_locations);

            % log(RRactual)
            %%%%%%%%%%%%%%%
            featMat_ldc = log(colvec(RR_intervals(this_iter_QRS_seq_idx))/obj.sampling_rate);

            % Update point
            obj.progress_hdl.checkpoint('Features calculation');

            % log(RRpost)
            %%%%%%%%%%%%%%%

            featMat_ldc = [ featMat_ldc log(colvec(RR_intervals(min(cant_QRS_locations, this_iter_QRS_seq_idx+1)))/obj.sampling_rate)];

            % Update point
            obj.progress_hdl.checkpoint('Features calculation');

            % log(RRmean1)
            %%%%%%%%%%%%%%%

            aux_idx = arrayfun(@(a)( find(  QRS_locations >= (QRS_locations(a) - 1*60*obj.sampling_rate) & ... 
                                            QRS_locations <= QRS_locations(a) )), ...
                               this_iter_seq_idx, 'UniformOutput', false);

            aux_featval = cellfun(@(a)(mean(RR_intervals(a))/obj.sampling_rate), aux_idx);
            featMat_ldc = [ featMat_ldc colvec(log(aux_featval)) ];

            % Update point
            obj.progress_hdl.checkpoint('Features calculation');

            % log(RRmean20)
            %%%%%%%%%%%%%%%

            aux_idx = arrayfun(@(a)( find(  QRS_locations >= (QRS_locations(a) - 20*60*obj.sampling_rate) & ... 
                                            QRS_locations <= QRS_locations(a) )), ...
                               this_iter_QRS_seq_idx, 'UniformOutput', false);

            aux_featval = cellfun(@(a)(mean(RR_intervals(a))/obj.sampling_rate), aux_idx);
            featMat_ldc = [ featMat_ldc colvec(log(aux_featval)) ];

            % Update point
            obj.progress_hdl.checkpoint('Features calculation');

            % AutoCorr1PlanoWTZeroCross
            %%%%%%%%%%%%%%%
            % AutoCorr1PlanoWTModMaxPos
            %%%%%%%%%%%%%%%

            aux_idx = arrayfun(@(a)( max(1, QRS_locations(a) - round(0.08*obj.sampling_rate)): ...
                                     min(this_iter_ECG_resampled_size, QRS_locations(a) + round(0.08*obj.sampling_rate))) , ...
                               this_iter_seq_idx, 'UniformOutput', false);

            aux_featval = cellfun(@(a)(squeeze(wtECG(a,:,obj.scale_idx(4)))), aux_idx, 'UniformOutput', false);

            %calculate PCA matrix in this slices for feature
            %AutoCorr1PlanoWTModMaxPos
            autovec_slices = cellfun(@(a)(autovec_calculation(a)), aux_featval, 'UniformOutput', false );

            aux_idx = arrayfun(@(a)( max(1, QRS_locations(a) - round(0.13*obj.sampling_rate)): ...
                                     min(this_iter_ECG_resampled_size, QRS_locations(a) + round(0.2*obj.sampling_rate))) , ...
                               this_iter_seq_idx, 'UniformOutput', false);

            aux_featval = cellfun(@(a,b)(squeeze(wtECG(a,:,obj.scale_idx(4))) * b), aux_idx, autovec_slices, 'UniformOutput', false);
            [ aux_mp aux_zc ] = cellfun(@(a)(CalcModMaxPos(a(:,1))), aux_featval );

            featMat_ldc = [ featMat_ldc colvec(aux_zc) colvec(aux_mp) ];

            % Update point
            obj.progress_hdl.checkpoint('Features calculation');

            % AutoCorr2PlanoWTZeroCross
            %%%%%%%%%%%%%%%
            % AutoCorr2PlanoWTModMaxPos
            %%%%%%%%%%%%%%%

            [ aux_mp aux_zc ] = cellfun(@(a)(CalcModMaxPos(a(:,2))), aux_featval );

            featMat_ldc = [ featMat_ldc colvec(aux_zc) colvec(aux_mp) ];


            % Update point
            obj.progress_hdl.checkpoint('Features calculation');

            %clear some auxiliar variables.
            clear aux*

            %calculate clustering features.


            % log(RRactual)
            %%%%%%%%%%%%%%%
            featMat_clust = featMat_ldc(:,1);

            % Update point
            obj.progress_hdl.checkpoint('Features calculation');

            % log(RRant)
            %%%%%%%%%%%%%%%
            featMat_clust = [ featMat_clust log(colvec(RR_intervals(max(1, this_iter_QRS_seq_idx-1)))/obj.sampling_rate)];

            % Update point
            obj.progress_hdl.checkpoint('Features calculation');

            % log(Prematuridad_loc)
            %%%%%%%%%%%%%%%

            aux_idx = arrayfun(@(a)(max(1,a-1):min(cant_QRS_locations,a+1)), this_iter_QRS_seq_idx, 'UniformOutput', false);
            %tengo que contemplar los casos extremos
            if( length(aux_idx{1}) < 3 )
                aux_idx{1} = [1 aux_idx{1}];
            end
            if( length(aux_idx{end}) < 3 )
                aux_idx{end} = [aux_idx{end} cant_QRS_locations];
            end

            aux_featval = cellfun(@(a)( RR_intervals(a(2))/sum(RR_intervals(a)) ), aux_idx);
            featMat_clust = [ featMat_clust colvec(log(aux_featval)) ];

            % Update point
            obj.progress_hdl.checkpoint('Features calculation');

            % log(AbsRRlocalVar)
            %%%%%%%%%%%%%%%

            aux_idx = arrayfun(@(a)(max(1,a-1):min(cant_QRS_locations,a+1)), this_iter_QRS_seq_idx, 'UniformOutput', false);
            if( length(aux_idx{1}) < 3 )
                aux_idx{1} = [1 aux_idx{1}];
            end
            if( length(aux_idx{end}) < 3 )
                aux_idx{end} = [aux_idx{end} cant_QRS_locations];
            end

            aux_featval = cellfun(@(a)(obj.SmallValue+sum(abs(dRR_intervals(a)))/obj.sampling_rate), aux_idx);

            featMat_clust = [ featMat_clust colvec(log(aux_featval)) ];

            % Update point
            obj.progress_hdl.checkpoint('Features calculation');

            % log(RRmean20)
            %%%%%%%%%%%%%%%

            featMat_clust = [ featMat_clust featMat_ldc(:,4) ];

            % Update point
            obj.progress_hdl.checkpoint('Features calculation');

            % log(QRS_max_scale_proj12)
            %%%%%%%%%%%%%%%

            aux_idx = arrayfun(@(a)( max(1, QRS_locations(a) - round(0.08*obj.sampling_rate)): ...
                                     min(this_iter_ECG_resampled_size, QRS_locations(a) + round(0.08*obj.sampling_rate))) , ...
                               this_iter_seq_idx, 'UniformOutput', false);

            aux_featval = cellfun(@(a)(wtECG(a,:,:)), aux_idx, 'UniformOutput', false);

            % Update point
            obj.progress_hdl.checkpoint('Features calculation');

            aux_featval = cellfun(@(a)(squeeze(sum(abs(a(:,1,:))))'), aux_featval, 'UniformOutput', false );
            aux_featval = cellfun(@(a)(obj.scales * colvec(a) ./ sum(a)), aux_featval );

        %     iAreaAbs = squeeze(sum(abs(wtAux)))';
        %     MaxProjArea = obj.scales * iAreaAbs ./ sum(iAreaAbs);

            featMat_clust = [ featMat_clust colvec(log(aux_featval)) ];

            % Update point
            obj.progress_hdl.checkpoint('Features calculation');

            % AutoCorr1PlanoWTModMaxPos
            %%%%%%%%%%%%%%%

            aux_idx = arrayfun(@(a)( max(1, QRS_locations(a) - round(0.13*obj.sampling_rate)): ...
                                     min(this_iter_ECG_resampled_size, QRS_locations(a) + round(0.2*obj.sampling_rate))) , ...
                               this_iter_seq_idx, 'UniformOutput', false);

            aux_featval = cellfun(@(a,b)(squeeze(wtECG(a,:,obj.scale_idx(4))) * b), aux_idx, autovec_slices, 'UniformOutput', false);
            aux_featval = cellfun(@(a)(CalcModMaxPos(a(:,1))), aux_featval);

            featMat_clust = [ featMat_clust colvec(aux_featval) ];


            % Update point
            obj.progress_hdl.checkpoint('Features calculation');

            % CrossCorr_QRST_ModMaxVal13
            %%%%%%%%%%%%%%%

            aux_idx = arrayfun(@(a)( max(1, QRS_locations(a) - round(0.13*obj.sampling_rate)): ...
                                     min( [ this_iter_ECG_resampled_size, ...
                                            QRS_locations(a) + round(0.4*obj.sampling_rate), ...
                                            QRS_locations(a) + RR_intervals(a) - round(0.13*obj.sampling_rate) ] ) ), ...
                               this_iter_seq_idx, 'UniformOutput', false);


            if( this_header.nsig > 2 ) 
                % wtECG already projected in obj.autovec.
                aux_featval = cellfun(@(a)( CalcModMaxVal( squeeze(wtECG(a,:,obj.scale_idx(3))) ) ), aux_idx);
            else
                aux_featval = cellfun(@(a)( CalcModMaxVal( squeeze(wtECG(a,:,obj.scale_idx(3))) * obj.autovec ) ), aux_idx);
            end

            featMat_clust = [ featMat_clust colvec(aux_featval) ];

            %save ECG for user interface and cluster analysis
            for ii = 1:length(obj.fnPayload)
                payload.(obj.fnPayload{ii}) = eval(obj.fnPayload{ii});
            end
    
        end
        
        function allPayloads = Concatenate(obj, allPayloads, payload )
            
            if( isempty(allPayloads) )
                for ii = 1:length(obj.fnFMfile)
                     allPayloads.(obj.fnFMfile{ii}) = payload.(obj.fnFMfile{ii});
                end
            else
                for ii = 1:length(obj.fnFMfile)
                     allPayloads.(obj.fnFMfile{ii}) = [allPayloads.(obj.fnFMfile{ii}); payload.(obj.fnFMfile{ii})];
                end
            end
            
            [~, rec_filename] = fileparts(obj.recording_name);

            for ii = 1:length(obj.fnECGfile)
                 aux_payload.(obj.fnECGfile{ii}) = payload.(obj.fnECGfile{ii});
            end
            
            save([obj.tmp_path 'tmpfile_' rec_filename '_ECG_cantpids_' num2str(obj.cant_pids) '_' num2str(obj.this_pid) '_' num2str(obj.this_iteration) '_' num2str(obj.cant_iteration)  '.mat'], '-struct', 'aux_payload' );

            aux_ECG_filename = [obj.tmp_path 'tmpfile_' rec_filename '_ECG_cantpids_' num2str(obj.cant_pids) '_thispid_' num2str(obj.this_pid) '_iteration_' num2str(obj.this_iteration) '_of_' num2str(obj.cant_iteration)  '.mat'];
            %do this for cluster supercomputers reasons
            movefile(   [obj.tmp_path 'tmpfile_' rec_filename '_ECG_cantpids_' num2str(obj.cant_pids) '_' num2str(obj.this_pid) '_' num2str(obj.this_iteration) '_' num2str(obj.cant_iteration)  '.mat'], ...
                aux_ECG_filename, 'f' );
            
            cant_QRS = length(payload.QRS_locations);
            obj.QRS_idx = [obj.QRS_idx; uint16(colvec(1:cant_QRS))];
            obj.File_idx = [obj.File_idx; repmat(uint16(obj.file_count), cant_QRS, 1)];
            obj.File_list = [obj.File_list; cellstr(aux_ECG_filename)];
            obj.file_count = obj.file_count + 1;

        end
        
        
        
    end
   
end
   
function autovec = autovec_calculation(wtECGslice)

    mean_wtECGslice = mean(wtECGslice);    
    wtECGslice_cov = cov( bsxfun( @minus, wtECGslice, mean_wtECGslice ));
    [autovec autoval] = eig(wtECGslice_cov); 
    [~, autoval_idx] = sort(diag(autoval), 'descend');
    autovec = autovec(:,autoval_idx);
    
end

function [ ModMaxPos ZeroCross ] = CalcModMaxPos( wtSignal )

    ModMaxPos = nan;
    ZeroCross = nan;

    lwtSignal = size(wtSignal, 1);

    if( all( wtSignal == 0 ) )
        return
    end

    ac1 = conv(wtSignal(:,1), flipud(wtSignal(:,1)));

    ac1 = ac1(lwtSignal:end);
    ac1 = 1/ac1(1)*ac1;

    interp_idx = (1:0.1:lwtSignal)';

    ac1_interp = spline( 1:lwtSignal, ac1, interp_idx );
    iZeroCrossPos_ac1 = myzerocros(ac1);

    if( isempty(iZeroCrossPos_ac1))
        return
    else
        if( iZeroCrossPos_ac1(1) < lwtSignal )
            if( sign(ac1(iZeroCrossPos_ac1(1))) ~= sign(ac1(iZeroCrossPos_ac1(1)+1)) )
                ZeroCross = iZeroCrossPos_ac1(1) - ac1(iZeroCrossPos_ac1(1)) / ( ac1(iZeroCrossPos_ac1(1)+1) - ac1(iZeroCrossPos_ac1(1)) );
            else
                ZeroCross = (iZeroCrossPos_ac1(1)-1) - ac1(iZeroCrossPos_ac1(1)-1) / ( ac1(iZeroCrossPos_ac1(1)) - ac1(iZeroCrossPos_ac1(1)-1) );
            end
        else
            return
        end   
    end

    ModMaxPos_ac1 = modmax( ac1_interp, 2, 0, 0 );

    if( ~isempty(ModMaxPos_ac1) )
        valid_ac1_max_idx = find(interp_idx(ModMaxPos_ac1) > ZeroCross);
        if( ~isempty(valid_ac1_max_idx) )
            ModMaxPos = interp_idx(ModMaxPos_ac1(valid_ac1_max_idx(1)));
        else
            return
        end
    else
        return
    end

end

function ModMaxVal = CalcModMaxVal( wtSignal )

    % std_wtSignal = std(wtSignal);
    % wtSignal = bsxfun( @times, wtSignal, 1./std_wtSignal);

    cc = conv(wtSignal(:,1), flipud(wtSignal(:,2)));

    [~, max_pos] = max(abs(cc));

    ModMaxVal = cc(max_pos);

end