classdef ECGtask_classification_features_calc < ECGtask

% ECGtask for ECGwrapper (for Matlab)
% ---------------------------------
% 
% Description:
% 
% Abstract class for defining ECGtask interface
% 
% 
% Author: Mariano Llamedo Soria (llamedom at frba.utn.edu.ar)
% Version: 0.1 beta
% Birthdate  : 18/2/2013
% Last update: 18/2/2013
       
    properties(GetAccess = public, Constant)
        name = 'ECG_hb_class_fc';
        target_units = 'ADCu';
        doPayload = true;
        % if user = memory;
        % memory_constant is the fraction respect to user.MaxPossibleArrayBytes
        % which determines the maximum input data size.

    end
    
    properties(GetAccess = public, SetAccess = private)
        true_labels
        memory_constant = 0.4;
        
        started = false;
        
    end
    
    properties( Access = private, Constant)
    
        sampling_rate = 360; %Hz
        fnPayload = {'featMat_clust' 'featMat_ldc' };
%         fnPayload = {'featMat_clust' 'featMat_ldc' 'ECG' 'QRS_locations' 'RR_intervals'};
        fnFMfile = {'featMat_clust' 'featMat_ldc'};
        SmallValue = 0.0001;
        
    end

    properties( Access = private)
        delayLP
        q_filters
        LP
        LPsize
        scale_idx
        scales    
        series
        cant_QRS_locations        
    end
    
    properties
        progress_handle 
        class_labeling
        autovec
        user_string
        tmp_path       

    end
    
    methods
           
        function obj = ECGclassificationTask(obj)
            % not implemented
        end
        
        function Start(obj, ECG_header, ECG_annotations)
            
            % the dimension of the input block should be affected for the
            % internal sampling rate, since a resampling is performed
            % internally.
            obj.memory_constant = obj.memory_constant * ECG_header.freq / 360;
            
            %% Calculate accesories signals
            
            %load low-pass preprocessing filter.
            load( ['LP Fs ' num2str(obj.sampling_rate) 'Hz - Fc 35Hz - 80db stop.mat' ]);
            obj.LPsize = length(LP.numerator);
            obj.delayLP = round((obj.LPsize-1)/2);
            obj.LP = LP;
            
            % Wavelet transform filter bank design
            obj.scales = 3:6;
            CantScales = length(obj.scales);
            MaxScales = max(obj.scales);
            obj.scale_idx= nan(MaxScales,1);

            filters_cache_filename = ['wt_filters_' num2str(obj.scales) ' scales_' num2str(obj.sampling_rate) ' Hz.mat' ];
            if( exist(filters_cache_filename, 'file') )
                aux = load( filters_cache_filename );
                obj.q_filters = aux.q_filters;
            else
                obj.q_filters = qs_filter_design(obj.scales, obj.sampling_rate);
            end
            
            for ii = 1:CantScales
                %indice para saber que escala corresponde a cada columna de MyWavelets
                obj.scale_idx(obj.scales(ii)) = ii;
            end

            if( isfield( ECG_annotations, 'anntyp'))
                obj.true_labels = ECG_annotations.anntyp;
            end
            QRS_locations = ECG_annotations.time;
            
            % RR interval sequence
            RR_intervals = diff(QRS_locations, 1);
            RR_intervals = colvec([ RR_intervals(1); RR_intervals ]);

            dRR_intervals = diff(RR_intervals, 1);
            dRR_intervals = colvec([ dRR_intervals(2); dRR_intervals(2); dRR_intervals(2:end) ]);

            % interval resampling.
            obj.series.RR_intervals = round(RR_intervals * obj.sampling_rate / ECG_header.freq);
            obj.series.dRR_intervals = round(dRR_intervals * obj.sampling_rate / ECG_header.freq);
            obj.series.QRS_locations = QRS_locations;

            obj.cant_QRS_locations = length(obj.series.RR_intervals);
            
            obj.started = true;
            
        end
        
        function payload = Process(obj, ECG, ECG_start_offset, ECG_sample_start_end_idx, ECG_header, ECG_annotations, ECG_annotations_start_end_idx )
            
            payload = [];
            
            %refered to the whole ECG recording
            this_iter_QRS_idx = ECG_annotations_start_end_idx(1):ECG_annotations_start_end_idx(2);
            RR_intervals = obj.series.RR_intervals(this_iter_QRS_idx);
            
            %refered to this ECG strip
            QRS_locations = ECG_annotations.time;
            
            %% Preprocessing
            
            % Update point
            obj.progress_handle.checkpoint('Resampling');
            
            %resample to obj.sampling_rate Hz
            wtECG = int16(round(resample(double(ECG), obj.sampling_rate, ECG_header.freq)));
            clear ECG
            QRS_locations = colvec(round(QRS_locations * obj.sampling_rate / ECG_header.freq));
        %         this_iter_true_labels = true_labels(this_iter_QRS_idx);
            this_iter_ECG_resampled_size = size(wtECG,1);

            if( ECG_header.nsig > 2 ) 
                %multilead approach.
                % project ECG and wtECG to obtain PCA components
                wtECG = int16(round(double(wtECG) * obj.autovec(:,1:2)));
        %             wtECG = cellfun(@(a)( a * obj.autovec(:,1:2) ), mat2cell(wtECG, this_iter_ECG_resampled_size, obj.ECG_header.nsig, ones(1,CantScales) ), 'UniformOutput', false);
        %             wtECG = cell2mat(wtECG);
            end

            % Update point
            obj.progress_handle.checkpoint('Filtering');

            %Low pass filtering @ 35Hz => ECG recovery
            wtECG = int16(round(filter(obj.LP, double(wtECG))));
            %delay compensation.
            wtECG = [ int16(zeros(obj.delayLP, size(wtECG,2) )) ;wtECG(obj.LPsize:end,:) ; int16(zeros(obj.delayLP, size(wtECG,2)))];

            %Quito la linea de base.
        %     ECG = BaselineWanderRemovalMedian( ECG, obj.sampling_rate);
            wtECG = int16(round(BaselineWanderRemovalSplines( double(wtECG), QRS_locations, obj.sampling_rate)));

            % Update point
            obj.progress_handle.checkpoint('Wavelet transform calculation');

            wtECG = int16(round(qs_wt(double(wtECG), obj.scales, obj.sampling_rate, obj.q_filters)));
            
            
            
            %% Features Calculation

            % Update point
            obj.progress_handle.checkpoint('RRi calculation');
            
            this_iter_cant_QRS = length(QRS_locations);
            this_iter_seq_idx = 1:this_iter_cant_QRS;

            % log(RRactual)
            %%%%%%%%%%%%%%%
            featMat_ldc = log(colvec(RR_intervals)/obj.sampling_rate);

            obj.progress_handle.checkpoint('RRi+1 calculation');

            % log(RRpost)
            %%%%%%%%%%%%%%%

            featMat_ldc = [ featMat_ldc log(colvec(obj.series.RR_intervals(min(obj.cant_QRS_locations, this_iter_QRS_idx+1)))/obj.sampling_rate)];

            % Update point
            obj.progress_handle.checkpoint('RRmean1 calculation');

            % log(RRmean1)
            %%%%%%%%%%%%%%%

            aux_idx = arrayfun(@(a)( findStartEnd(  obj.series.QRS_locations >= (obj.series.QRS_locations(a) - 1*60*obj.sampling_rate) & ... 
                                                    obj.series.QRS_locations <= obj.series.QRS_locations(a) )), ...
                               this_iter_QRS_idx, 'UniformOutput', false);

            aux_featval = cellfun(@(a)(mean(obj.series.RR_intervals(a(1):a(2)))/obj.sampling_rate), aux_idx);
            featMat_ldc = [ featMat_ldc colvec(log(aux_featval)) ];

            % Update point
            obj.progress_handle.checkpoint('RRmean20 calculation');

            % log(RRmean20)
            %%%%%%%%%%%%%%%

            aux_idx = arrayfun(@(a)( findStartEnd(  obj.series.RR_intervals >= (obj.series.RR_intervals(a) - 20*60*obj.sampling_rate) & ... 
                                                    obj.series.RR_intervals <= obj.series.RR_intervals(a) )), ...
                               this_iter_QRS_idx, 'UniformOutput', false);

            aux_featval = cellfun(@(a)(mean(obj.series.RR_intervals(a(1):a(2)))/obj.sampling_rate), aux_idx);
            featMat_ldc = [ featMat_ldc colvec(log(aux_featval)) ];

            % Update point
            obj.progress_handle.checkpoint('AC1 calculation');

            % AutoCorr1PlanoWTZeroCross
            %%%%%%%%%%%%%%%
            % AutoCorr1PlanoWTModMaxPos
            %%%%%%%%%%%%%%%

            aux_idx = arrayfun(@(a)( [ max(1, QRS_locations(a) - round(0.08*obj.sampling_rate))  ...
                                     min(this_iter_ECG_resampled_size, QRS_locations(a) + round(0.08*obj.sampling_rate)) ]) , ...
                               this_iter_seq_idx , 'UniformOutput', false);

            %calculate PCA matrix in this slices for feature
            %AutoCorr1PlanoWTModMaxPos

            aux_idx2 = arrayfun(@(a)( [ max(1, QRS_locations(a) - round(0.13*obj.sampling_rate))  ...
                                     min(this_iter_ECG_resampled_size, QRS_locations(a) + round(0.2*obj.sampling_rate)) ]) , ...
                               this_iter_seq_idx, 'UniformOutput', false);

            aux_featval = cellfun(@(a,b)(double(squeeze(wtECG(b(1):b(2),:,obj.scale_idx(4)))) * autovec_calculation(double(squeeze(wtECG(a(1):a(2),:,obj.scale_idx(4))))) ), aux_idx, aux_idx2, 'UniformOutput', false);
            
            [ aux_mp aux_zc ] = cellfun(@(a)(CalcModMaxPos(a(:,1))), aux_featval );

            featMat_ldc = [ featMat_ldc colvec(aux_zc) colvec(aux_mp) ];

            % Update point
            obj.progress_handle.checkpoint('AC2 calculation');

            % AutoCorr2PlanoWTZeroCross
            %%%%%%%%%%%%%%%
            % AutoCorr2PlanoWTModMaxPos
            %%%%%%%%%%%%%%%

            [ aux_mp aux_zc ] = cellfun(@(a)(CalcModMaxPos(a(:,2))), aux_featval );

            featMat_ldc = [ featMat_ldc colvec(aux_zc) colvec(aux_mp) ];

            % Update point
            obj.progress_handle.checkpoint('RRi replication');

            %calculate clustering features.


            % log(RRactual)
            %%%%%%%%%%%%%%%
            featMat_clust = featMat_ldc(:,1);

            % Update point
            obj.progress_handle.checkpoint('RRi-1 calculation');

            % log(RRant)
            %%%%%%%%%%%%%%%
            featMat_clust = [ featMat_clust log(colvec(obj.series.RR_intervals(max(1, this_iter_QRS_idx-1)))/obj.sampling_rate)];

            % Update point
            obj.progress_handle.checkpoint('RRp calculation');

            % log(Prematuridad_loc)
            %%%%%%%%%%%%%%%

            aux_idx = arrayfun(@(a)(max(1,a-1):min(obj.cant_QRS_locations,a+1)), this_iter_QRS_idx, 'UniformOutput', false);
            %tengo que contemplar los casos extremos
            if( length(aux_idx{1}) < 3 )
                aux_idx{1} = [1 aux_idx{1}];
            end
            if( length(aux_idx{end}) < 3 )
                aux_idx{end} = [aux_idx{end} obj.cant_QRS_locations];
            end

            aux_featval = cellfun(@(a)( obj.series.RR_intervals(a(2))/sum(obj.series.RR_intervals(a)) ), aux_idx);
            featMat_clust = [ featMat_clust colvec(log(aux_featval)) ];

            % Update point
            obj.progress_handle.checkpoint('RRlv calculation');

            % log(AbsRRlocalVar)
            %%%%%%%%%%%%%%%

            aux_idx = arrayfun(@(a)(max(1,a-1):min(obj.cant_QRS_locations,a+1)), this_iter_QRS_idx, 'UniformOutput', false);
            if( length(aux_idx{1}) < 3 )
                aux_idx{1} = [1 aux_idx{1}];
            end
            if( length(aux_idx{end}) < 3 )
                aux_idx{end} = [aux_idx{end} obj.cant_QRS_locations];
            end

            aux_featval = cellfun(@(a)(obj.SmallValue+sum(abs(obj.series.dRR_intervals(a)))/obj.sampling_rate), aux_idx);

            featMat_clust = [ featMat_clust colvec(log(aux_featval)) ];

            % Update point
            obj.progress_handle.checkpoint('RRmean20 replication');

            % log(RRmean20)
            %%%%%%%%%%%%%%%

            featMat_clust = [ featMat_clust featMat_ldc(:,4) ];

            % Update point
            obj.progress_handle.checkpoint('QRSscale calculation');

            % log(QRS_max_scale_proj12)
            %%%%%%%%%%%%%%%

            aux_idx = arrayfun(@(a)( [ max(1, QRS_locations(a) - round(0.08*obj.sampling_rate)) ...
                                     min(this_iter_ECG_resampled_size, QRS_locations(a) + round(0.08*obj.sampling_rate)) ]) , ...
                               this_iter_seq_idx, 'UniformOutput', false);

            aux_featval = cellfun(@(a)( Calc_QRS_max_scale_proj( double(wtECG(a(1):a(2),:,:)), obj.scales) ), aux_idx);

        %     iAreaAbs = squeeze(sum(abs(wtAux)))';
        %     MaxProjArea = scales * iAreaAbs ./ sum(iAreaAbs);

            featMat_clust = [ featMat_clust colvec(log(aux_featval)) ];

            % Update point
            obj.progress_handle.checkpoint('AC1 replication');

            % AutoCorr1PlanoWTModMaxPos
            %%%%%%%%%%%%%%%

            featMat_clust = [ featMat_clust featMat_ldc(:,6) ];

            % Update point
            obj.progress_handle.checkpoint('CC calculation');

            % CrossCorr_QRST_ModMaxVal13
            %%%%%%%%%%%%%%%

            aux_idx = arrayfun(@(a)( [max(1, QRS_locations(a) - round(0.13*obj.sampling_rate)) ...
                                     min( [ this_iter_ECG_resampled_size, ...
                                            QRS_locations(a) + round(0.4*obj.sampling_rate), ...
                                            QRS_locations(a) + RR_intervals(a) - round(0.13*obj.sampling_rate) ] ) ]), ...
                               this_iter_seq_idx, 'UniformOutput', false);

            if( ECG_header.nsig > 2 ) 
                % wtECG already projected in autovec.
                aux_featval = cellfun(@(a)( CalcModMaxVal( double(squeeze(wtECG(a(1):a(2),:,obj.scale_idx(3)))) ) ), aux_idx);
            else
                aux_featval = cellfun(@(a)( CalcModMaxVal( double(squeeze(wtECG(a(1):a(2),:,obj.scale_idx(3)))) * obj.autovec ) ), aux_idx);
            end
            
            featMat_clust = [ featMat_clust colvec(aux_featval) ];

            %save ECG for user interface and cluster analysis
            for ii = 1:length(obj.fnPayload)
                payload.(obj.fnPayload{ii}) = eval(obj.fnPayload{ii});
            end

        end
        
        function payload = Finish(obj, payload, ECG_header)
            % not implemented
        end
        
        function plA = Concatenate(obj, plA, plB)
            
            if( isempty(plA) )
                for ii = 1:length(obj.fnFMfile)
                     plA.(obj.fnFMfile{ii}) = plB.(obj.fnFMfile{ii});
                end
            else
                for ii = 1:length(obj.fnFMfile)
                     plA.(obj.fnFMfile{ii}) = [plA.(obj.fnFMfile{ii}); plB.(obj.fnFMfile{ii})];
                end
            end
                        
        end

    end
    
    methods ( Access = private )
        
        
    end
    
end


%%%%%%%%%%%%%%%%%%
% Other functions
%%%%%%%%%%%%%%%%%%


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

function start_end_aux = findStartEnd( bAux )

    start_aux = find(  bAux, 1, 'first' );
    end_aux = find(  bAux, 1, 'last' );
    start_end_aux = [start_aux end_aux];
                                        
end

function aux_slice = Calc_QRS_max_scale_proj(aux_slice, scales)
    aux_slice = squeeze(sum(abs(aux_slice(:,1,:))))';
    aux_slice = scales * colvec(aux_slice) ./ sum(aux_slice);
end
