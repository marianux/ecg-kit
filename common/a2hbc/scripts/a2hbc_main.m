function [ Labels, ConfusionMatrix, LabelList ] = a2hbc_main(varargin)

% A2HBC main script
% -----------------
% See a2hbc.m in the parent directory or the documentation for help.
% 
% Author: Mariano Llamedo Soria (llamedom at frba.utn.edu.ar)
% Birthdate  : 16/8/2011
% Last update: 9/2/2012

%% Constants and definitions
persistent userNoAnswers

if( isempty(userNoAnswers) )
    %first use
    userNoAnswers = 1;
end

a2hbc_header;

% some variable declarations;
Labels = [];
ConfusionMatrix = [];
LabelList = [];
true_labels = [];
bArgCheckPassed = false;
CantSamplesClose = 10;
CantSamplesFar = 10;
bUserExit = false;
bRecordingChanged = true;
bLabelingChanged = false;
repeat_idx = 1;

%% Argument parsing

%argument definition
p = inputParser;   % Create instance of inputParser class.
p.addParamValue('recording_name', [], @(x)( ischar(x)));
p.addParamValue('recording_format', [], @(x)( ischar(x) && any(strcmpi(x,cKnownFormats))) );
p.addParamValue('ECG', [], @(x)(isnumeric(x)) );
p.addParamValue('ECG_header', [], @(x)(isstruct(x)) );
p.addParamValue('ECG_annotations', [], @(x)( isstruct(x) ) );
p.addParamValue('op_mode', 'auto', @(x)( (isnumeric(x) && x >= 1  && x <= maxOperMode) || any(strcmpi(x,cKnownModesOfOperation))) );
p.addParamValue('cant_pids', 1, @(x)(isnumeric(x) && x > 0 ) );
p.addParamValue('this_pid', 1, @(x)(isnumeric(x) && x > 0 ) );
p.addParamValue('CacheData', true, @(x)(islogical(x)) );
p.addParamValue('InteractiveMode', false, @(x)(islogical(x)) );
p.addParamValue('SimulateExpert', false, @(x)(islogical(x)) );
p.addParamValue('tmp_path', [], @(x)(ischar(x)) );
p.addParamValue('NumOfClusters', 12, @(x)(isnumeric(x) && x > 1 ) );
p.addParamValue('ClusteringRepetitions', 1, @(x)(isnumeric(x) && x > 0 && x <= 10 ) );
p.addParamValue('ClusterPresence', 75, @(x)(isnumeric(x) && x >= 0 && x <= 100 ) );
p.addParamValue('Repetitions', 1, @(x)(isnumeric(x) && x > 0 ) );
%additions
p.addParamValue('class_labeling', 'AAMI', @(x)( ischar(x) && any(strcmpi(x,cKnownLabelings))) );
p.addParamValue('class_filter', char({'Normal', 'Supraventricular', 'Ventricular'}), @(x)( ischar(x) || iscell(x) ) );


try
    p.parse( varargin{:} );
catch MyError
    rethrow(MyError);    
end

recording_name = p.Results.recording_name;
recording_format = p.Results.recording_format;
ECG_total = p.Results.ECG;
ECG_header = p.Results.ECG_header;
ECG_annotations = p.Results.ECG_annotations;
bCache = p.Results.CacheData;
op_mode = p.Results.op_mode;
cant_pids = p.Results.cant_pids;
this_pid = p.Results.this_pid;
tmp_path = p.Results.tmp_path;
bInteractiveMode = p.Results.InteractiveMode;
bSimulateExpert = p.Results.SimulateExpert;
CantClusters = p.Results.NumOfClusters;
iter_times = p.Results.ClusteringRepetitions;
cluster_presence = p.Results.ClusterPresence;
Repetitions = p.Results.Repetitions;
class_labeling = p.Results.class_labeling;
class_filter = p.Results.class_filter;

% Dont know why this variable uses a lot of bytes to store at disk.
clear p

%% Start of the algorithm

if( bHaveUserInterface )
    if( isempty(varargin) )
        bInteractiveMode = true;
    end
else
    bInteractiveMode = false;
end

typical_lablist = char(typical_lablists{find(strcmpi(class_labeling, cKnownLabelings))});
typical_cant_labels = size(typical_lablist,1);
this_iter_cm = nan(typical_cant_labels);
ConfusionMatrix = zeros(typical_cant_labels, typical_cant_labels,Repetitions);

while ( ~bUserExit && repeat_idx <= Repetitions )

    try

        %% work starts here
        
        if( bRecordingChanged )
            % new recording, needs processing
            bRecordingChanged = false;
            
            %% Argument check

            % op_mode parsing
            if( isnumeric(op_mode) )
                op_mode = cKnownModesOfOperation{op_mode};
            end

            % class filtering parsing
            class_filter = char(intersect( cellstr(class_filter), cellstr(typical_lablist)));
            if( isempty(class_filter) )
                strAux = [ repmat(' + ', size(typical_lablist,1), 1) char(typical_lablist) repmat('\n', size(typical_lablist,1), 1 ) ];
                error( 'a2hbc:ArgCheck:InvalidClassFilter', ['Invalid class-filter. Please provide one of these classes:\n' rowvec(strAux')] );
            end

            bECG_sample_provided = false;
            
            %ECG parsing
            if( ~isempty(ECG_total) )
                % ECG already read
                
                bECG_sample_provided = true;

                
            elseif( ~isempty(recording_name) )

                [~, rec_filename] = fileparts(recording_name);

                % ECG to be read, create a wrapper object, later assign tasks to it.
                ECGw = ECGwrapper( ... 
                                    'recording_name', recording_name, ...
                                    'recording_format', recording_format, ...
                                    'tmp_path', tmp_path, ...
                                    'this_pid', sprintf('%d/%d', this_pid, cant_pids) ...
                                    );
                                
                ECG_header = ECGw.ECG_header;
                ECG_annotations = ECGw.ECG_annotations;
                
                % annotations are already AAMI, this is in order not to
                % convert again.
                recording_format = '';
            else
            %     strAux = help('a2hbc'); %#ok<MCHLP>
                error( 'a2hbc:ArgCheck:InvalidECGarg', 'Please provide an ECG recording as described in the documentation, help(''a2hbc'') maybe could help you.\n' );
            end
            
            if( isempty(ECG_header))
               error( 'a2hbc:ArgCheck:InvalidHeader', 'Please provide the ECG header.\n\n' );
            else
                if( ~isfield(ECG_header, cHeaderFieldNamesRequired ) )
                    strAux = [ repmat(' + ', length(cHeaderFieldNamesRequired), 1) char(cHeaderFieldNamesRequired) repmat('\n', length(cHeaderFieldNamesRequired), 1 ) ];
                    error( 'a2hbc:ArgCheck:InvalidHeader', ['Please provide the following fields in the header struct:\n ' rowvec(strAux') ] );
                end
            end

            if( isempty(ECG_annotations))
               error( 'a2hbc:ArgCheck:InvalidAnnotations', 'Please provide the ECG annotations.\n\n' );
            else
                if( ~isfield(ECG_annotations, cAnnotationsFieldNamesRequired ) )
                    strAux = [ repmat(' + ', length(cAnnotationsFieldNamesRequired), 1) char(cAnnotationsFieldNamesRequired) repmat('\n', length(cAnnotationsFieldNamesRequired), 1 ) ];
                    error( 'a2hbc:ArgCheck:InvalidAnnotations', ['Please provide the following fields in the annotations struct:\n ' rowvec(strAux') ] );
                end
            end

            [QRS_locations, aux_val ] = Annotation_process(ECG_annotations, recording_format, class_labeling);

            if( isempty(QRS_locations))
               error( 'a2hbc:ArgCheck:InvalidAnnotations', 'No QRS annotations available, please provide valid annotations or perform QRS detection.\n\n' );
            end
            
            true_labels = nan(size(aux_val));
            for ii = 1:length(typical_lablists_anntyp)
                bAux = aux_val == typical_lablists_anntyp{ii};
                true_labels(bAux) = ii;
            end

            rec_filename = ECG_header.recname;
                
            bLabelingChanged = false;

            lablist_idx = find(strcmpi(class_labeling, cKnownLabelings));

            %Update the automatic classifier.
            global_classifier = global_classifiers{lablist_idx};

            if( isempty(tmp_path) )
                tmp_path = tempdir;
            end

            %check path integrity.
            if(~exist(tmp_path, 'dir'))
                %try to create it
                if( ~mkdir(tmp_path) )
                    error('a2hbc:ArgCheck:InvalidPath', 'Invalid tmp_path. Please provide a valid path.\n' );
                end
            end

            bArgCheckPassed = true;

            

            %% Calculate accesories signals and coefficients
            
            if( this_pid > cant_pids )
                %% Classification only PIDs
                % Wait other PIDS sync here after Master build the featmat
                % file.

                %last pid
                bContinue = true;
                
                CachedFeatMatFiles = dir([tmp_path 'tmpfile_a2hbc_' ECGtask_class_fc_hdl.name '_' rec_filename '*.mat']);
                
                %Wait for Time2WaitPIDs seconds the finalization of all PIDs. Otherwise exit
                %with error.
                Start2Wait = tic();

                while(bContinue)

                    try

                        if( isempty( CachedFeatMatFiles ) )
                            error('a2hbc:PIDnotFinished', 'Handled error');
                        else                            
                            % exit the waiting loop.
                            bContinue = false;
                        end

                    catch ME


                        if( strfind(ME.identifier, 'a2hbc') )
                            if( toc(Start2Wait) > 2*Time2WaitPIDs )
                                error('a2hbc:PIDnotFinished', 'Timeout. Classification-only Slave give up waitng Master.');
                            end
                            pause(60);
                        else
                            rethrow(ME)
                        end
                    end


                end
                
            end

            % one task for calculating the PCA basis over the entire
            % recording
            ECGtask_PCA_proj_basis_hdl = ECGtask_PCA_proj_basis();

            % new task for calculating the features
            ECGtask_class_fc_hdl = ECGtask_classification_features_calc();
            
            %% ECG processing
            
            CachedFeatMatFiles = dir([tmp_path 'tmpfile_a2hbc_' ECGtask_class_fc_hdl.name '_' rec_filename '*.mat']);

            if( isempty( CachedFeatMatFiles ) )
                
                if( bECG_sample_provided )

                    aux_struct.time = QRS_locations;
                    aux_struct.anntyp = true_labels;
                    
                    if( ECG_header.nsig > 1  )
                        % As no wrapper is needed, I manually execute the task
                        ECGtask_PCA_proj_basis_hdl.Start(ECG_header, aux_struct);
                        ECGtask_PCA_proj_basis_hdl.Process(ECG_total, 1, [1 ECG_header.nsig], ECG_header, aux_struct, [1 length(QRS_locations) ] );
                        ECGtask_PCA_proj_basis_hdl.Finish([], ECG_header);
                    else
                        % PCA not need.
                        ECGtask_PCA_proj_basis_hdl.autovec = 1;
                    end

                    % Activate the progress_struct bar.
                    pb = progress_bar(rec_filename);
                    
                    ECGtask_class_fc_hdl.progress_handle = pb;
                    ECGtask_class_fc_hdl.autovec = ECGtask_PCA_proj_basis_hdl.autovec;
                    ECGtask_class_fc_hdl.Start(ECG_header, aux_struct);
                    payload = ECGtask_class_fc_hdl.Process(ECG_total, 1, [1 ECG_header.nsamp], ECG_header, aux_struct, [1 length(QRS_locations) ] ); %#ok<NASGU>

                    CachedFeatMatFiles = [tmp_path 'tmpfile_a2hbc_' ECGtask_class_fc_hdl.name '_' rec_filename '.mat'];
                    save(CachedFeatMatFiles, '-struct', 'payload');
                    clear payload
                    CachedFeatMatFiles = dir(CachedFeatMatFiles);
                    
                    delete(pb)
                    
                else

                    aux_struct.time = QRS_locations;
                    aux_struct.anntyp = true_labels;
                    ECGw.ECG_annotations = aux_struct;
                    
                    if( ECGw.ECG_header.nsig > 1  )
    %                     if( isempty(ECGw.ECG_annotations) )
    %                         ECGtask_PCA_proj_basis_hdl.cant_QRS_locations = 0;
    %                     else
    %                         ECGtask_PCA_proj_basis_hdl.cant_QRS_locations = length(ECGw.ECG_annotations.time);
    %                     end

                        ECGw.partition_mode =  'QRS';
                        ECGw.ECGtaskHandle = ECGtask_PCA_proj_basis_hdl;

                        %calculate the basis
                        ECGw.Run();
                    else
                        % PCA not need.
                        ECGtask_PCA_proj_basis_hdl.autovec = 1;
                    end

                    ECGw.class_labeling = class_labeling;
                    ECGtask_class_fc_hdl.autovec = ECGtask_PCA_proj_basis_hdl.autovec;

                    % assign the task to the wrapper.
                    ECGw.ECGtaskHandle = ECGtask_class_fc_hdl;
                    ECGw.partition_mode =  'QRS';

                    %% Iterations over the whole ECG recording
                    ECGw.Run;

                    CachedFeatMatFiles = [];
                    [tmp_path, CachedFeatMatFiles.name] = fileparts(ECGw.Result_files{1});
                    tmp_path = [tmp_path filesep];
                end
                
            end
            

            if( isempty( CachedFeatMatFiles ) )
                error( 'a2hbc:FeaturesNotCalc', 'Cached features files not found. Check files in the temporary dir.\n' );
            else
                
                % Activate the progress_struct bar.
                pb = progress_bar(rec_filename);
                
                % Update point
                pb.checkpoint('Restoring cached feature matrix');
                
                if( length(CachedFeatMatFiles) > 1 )
                    warning('TODO: generalizar para featmatrices arb grandes. Leyendo solo una parte de la feat matrix');
                    CachedFeatMatFileName = CachedFeatMatFiles(1).name;
                else
                    CachedFeatMatFileName = CachedFeatMatFiles.name;
                end
                
                % restore cached data.
                load( [ tmp_path CachedFeatMatFileName]);
                
            end

            if( bECG_sample_provided )
                ann = ECG_annotations;
            else
                ann = ECGw.ECG_annotations;
                aux_val = ann.anntyp;
                
                true_labels = nan(size(aux_val));
                for ii = 1:length(typical_lablists_anntyp)
                    bAux = aux_val == typical_lablists_anntyp{ii};
                    true_labels(bAux) = ii;
                end
                
            end
            
            % clear unused objects before continue
%             clear ECGw
%             clear ECGtask*
            
%           recorrer ECGw.Result_files y armar las matrices de features.

            %% Feature Matrices building


            %clean featureMatrix from NAN values.
%             AnyNaN_idx = find(any(isnan(featMat_clust),2));
%             iAllMedian = nanmedian(featMat_clust);
%             if(any(isnan(iAllMedian)))
%                error( 'a2hbc:AllNanFeatures', 'Features not calculated in this recording. Check recording or ask help.\n' );
%             end
%             for ii = rowvec(AnyNaN_idx)
%                bIdx2 = find(isnan(featMat_clust(ii,:)));
%                featMat_clust(ii, bIdx2) = iAllMedian(bIdx2);
%             end
% 
%             AnyNaN_idx = find(any(isnan(featMat_ldc),2));
%             iAllMedian = nanmedian(featMat_ldc);
%             if(any(isnan(iAllMedian)))
%                error( 'a2hbc:AllNanFeatures', 'Features not calculated in this recording. Check recording or ask help.\n' );
%             end
%             for ii = rowvec(AnyNaN_idx)
%                bIdx2 = find(isnan(featMat_ldc(ii,:)));
%                featMat_ldc(ii, bIdx2) = iAllMedian(bIdx2);
%             end

            % por el momento sigo usando el PRtools.
            % armo datasets.
            featMat_clust = prdataset(featMat_clust);
            featMat_ldc = prdataset(featMat_ldc);

            % Update point
            %clean featureMatrix from NAN values.
            pb.checkpoint('Clean featureMatrix from NAN values');
            
            featMat_ldc = deNaN_dataset(featMat_ldc, 'change');
            featMat_clust = deNaN_dataset(featMat_clust, 'change');
            
            featMat_ldc = deNaN_dataset(featMat_ldc, 'change', @isinf);
            featMat_clust = deNaN_dataset(featMat_clust, 'change', @isinf);
            

            if( bHaveUserInterface && (~exist('UCP_struct', 'var') || ~ishandle(UCP_struct.fig_hdl)) )
                UserControlPanel;
                ParseUserControlPanelInput;
            end

        else
            %% Only labeling changed ?
            
            if( bLabelingChanged )
                
                bLabelingChanged = false;
                
                % labeling changed
                [QRS_locations, aux_val ] = Annotation_process(ann, recording_format, class_labeling);

                true_labels = nan(size(aux_val));
                for ii = 1:length(typical_lablists_anntyp)
                    bAux = aux_val == typical_lablists_anntyp{ii};
                    true_labels(bAux) = ii;
                end
                
                %change lablists
                lablist_idx = find(strcmpi(class_labeling, cKnownLabelings));
                typical_lablist = char(typical_lablists{lablist_idx});
                typical_cant_labels = size(typical_lablist,1);
                this_iter_cm = nan(typical_cant_labels);
                ConfusionMatrix = zeros(typical_cant_labels, typical_cant_labels,Repetitions);
                
                %Update the automatic classifier.
                global_classifier = global_classifiers{lablist_idx};
                
            end
            
        end

        
        %% Classification

        cant_QRS_locations = size( featMat_ldc,1);
        
        already_labeled_idx = [];
        pending_hb_idx = setdiff((1:cant_QRS_locations)', already_labeled_idx );
        bCancel = false;

        % nasty addons
        RR_intervals = diff(QRS_locations, 1);
        RR_intervals = colvec([ RR_intervals(1); RR_intervals ]);
        
        if( ~bECG_sample_provided )
            ECG_total = ECGw.read_signal(1, ECGw.ECG_header.nsamp);
        end
        
        if( bHaveUserInterface )
            DisplayConfiguration;
        end
        
        while( ~bCancel && ~isempty(pending_hb_idx) )

            bDataSkiped = false;
            this_loop_already_labeled_idx = [];
            this_loop_cant_pending = length(pending_hb_idx);

            % Update point
            pb.checkpoint('Clustering data.');

            %skip warnings here
            warning off all; prwarning(0);
            Clust_Labels = cluster_data_with_EM_clust(featMat_clust(pending_hb_idx,:), qdc_new([],1e-6,1e-6, []), CantClusters, iter_times);
            warning on all; prwarning(1);

            [~, sort_idx] = sort( Clust_Labels );
            [all_clusters, aux_location] = unique(Clust_Labels(sort_idx,:), 'first');

            aux_location = [colvec(aux_location); cant_QRS_locations+1];
            cluster_sizes = diff(aux_location);
            cant_clusters = size(all_clusters,1);

            if( bHaveUserInterface)
                set(UCP_struct.Axes_hdl, 'Visible','on' );
%                 waitbar(0, UCP_struct.fig_hdl, 'Classification progress 0%')
                hb_percent_done = (length(already_labeled_idx) + length(this_loop_already_labeled_idx))/cant_QRS_locations;
                waitbar( hb_percent_done, UCP_struct.fig_hdl, [ 'Classification progress ' num2str(round(hb_percent_done*100)) '%' ])
            end
            
            %% Intento de clasificacion autom�tica.
            % Aqu� trabaja el GC.

            % Update point
            pb.checkpoint('Automatic classification');

            if( strcmpi(op_mode, 'slightly-assisted') || strcmpi(op_mode, 'auto') )

                dsThisPatient_classification = featMat_ldc * global_classifier;

                dsThisPatient_classification = labeld(dsThisPatient_classification);
                dsThisPatient_classification = renumlab(dsThisPatient_classification, typical_lablist); % try to fit on original lablist

                for ii = 1:cant_clusters

                    %los latidos de este cluster
                    aux_idx = sort_idx(aux_location(ii):(aux_location(ii+1)-1));

                    %busco las clasificaciones no rechazadas.
                    aux_idx1 = find(dsThisPatient_classification(aux_idx) ~= 0);

                    cantXclase = histc( dsThisPatient_classification(aux_idx(aux_idx1)), 1:typical_cant_labels );

                    %busco la clase mas predecida
                    [~, class_idx] = max(cantXclase);

                    if( cantXclase(class_idx) > ( cluster_presence  * 0.01 * cluster_sizes(ii) ) )
                        % La mayoria impone su clase.
                        Clust_Labels(aux_idx) = class_idx;
                        
                        this_loop_already_labeled_idx = [ this_loop_already_labeled_idx; colvec(aux_idx) ];

                        if( bHaveUserInterface)                            
                            hb_percent_done = (length(already_labeled_idx) + length(this_loop_already_labeled_idx))/cant_QRS_locations;
                            waitbar( hb_percent_done, UCP_struct.fig_hdl, [ 'Classification progress ' num2str(round(hb_percent_done*100)) '%' ])
                        end
                        
                    else
                        %for the slightly assisted, in the next section the
                        %expert will do the labeling.
                        if( strcmpi(op_mode, 'auto') )
                            %se deja la clasificacion autom�tica.
                            Clust_Labels(aux_idx) = dsThisPatient_classification(aux_idx);
                            
                            this_loop_already_labeled_idx = [ this_loop_already_labeled_idx; colvec(aux_idx) ];

                            if( bHaveUserInterface)                            
                                hb_percent_done = (length(already_labeled_idx) + length(this_loop_already_labeled_idx))/cant_QRS_locations;
                                waitbar( hb_percent_done, UCP_struct.fig_hdl, [ 'Classification progress ' num2str(round(hb_percent_done*100)) '%' ])
                            end                                
                            
                        end
                    end
                end                    
            end

            %% Clasificacion manual.
            % Aca consultamos al ORACLE.

            % Update point
            pb.checkpoint('Expert assistance');

            if( any(strcmpi(op_mode, {'assisted' 'slightly-assisted'}) ) )
                %Clustering + global classifier (GC) + oracle
                %Clustering + oracle
                %busco los centroides o algun elemento
                %representativo de cada cluster
                not_labeled_idx = setdiff(1:this_loop_cant_pending, this_loop_already_labeled_idx );
                [~, sort_idx] = sort( Clust_Labels(not_labeled_idx) );
                [all_clusters, aux_location] = unique(Clust_Labels(not_labeled_idx(sort_idx)), 'first');

                m_not_labeled = length(not_labeled_idx);

                aux_location = [colvec(aux_location); m_not_labeled+1];
                cluster_sizes = diff(aux_location);
                cant_clusters = size(all_clusters,1);

                ii = 1;
                while( ii <= cant_clusters )

                    %los latidos de este cluster
                    aux_idx = sort_idx(aux_location(ii):(aux_location(ii+1)-1));

                    %si es muy grande lo sampleo uniformemente para
                    %no generar mucho costo computacional.
                    laux_idx = length(aux_idx);
                    if( laux_idx > 5000 )
                        aux_idx2 = randsample(laux_idx, 5000);
                    else
                        aux_idx2 = 1:laux_idx;
                    end

                    %busco el centroide del cluster para
                    %etuiquetarlo.
                    clust_data_m = length(aux_idx2);
                    
                    if( clust_data_m > 1 )
                        
                        clust_data = [ featMat_clust(pending_hb_idx(not_labeled_idx(aux_idx(aux_idx2))),:); featMat_clust(pending_hb_idx(not_labeled_idx(aux_idx(aux_idx2))),:) ];
                        clust_data = setlabels(clust_data,[ones(clust_data_m,1);2*ones(clust_data_m,1)] );
                        clust_data = setprior(clust_data,[1 1]./2 );
                        clust_data_model = qdc(clust_data, 1e-6, 1e-6);
                        clust_data_posterior = clust_data * clust_data_model;
                        clust_dist2mu = +clust_data_posterior(1:clust_data_m,1);

                        [~, clust_dist2mu_sorted_idx ] = sort(clust_dist2mu, 'descend');

                    else
                        clust_dist2mu_sorted_idx = 1;
                    end
                    
                    if( ~bHaveUserInterface || bSimulateExpert )

                        clust_sample_centroid_idx = clust_dist2mu_sorted_idx(1);
                        %a cada cluster le asignamos la etiqueta
                        %VERDADERA del elemento mas cercano al centroide.
                        Clust_Labels(not_labeled_idx(aux_idx)) = true_labels(not_labeled_idx(aux_idx(aux_idx2(clust_sample_centroid_idx))));
                        
                        this_loop_already_labeled_idx = [ this_loop_already_labeled_idx; colvec(not_labeled_idx(aux_idx)) ];

                        if( bHaveUserInterface )
                            hb_percent_done = (length(already_labeled_idx) + length(this_loop_already_labeled_idx))/cant_QRS_locations;
                            waitbar( hb_percent_done, UCP_struct.fig_hdl, [ 'Classification progress ' num2str(round(hb_percent_done*100)) '%' ])
                        end
                        
                        ii = ii + 1;
                        
                    else
                        clust_sample_centroid_idx = clust_dist2mu_sorted_idx(1);
                        fprintf(1, ['Suggestion: ' typical_lablist(true_labels(not_labeled_idx(aux_idx(aux_idx2(clust_sample_centroid_idx)))),:) '\n']);
                        
                        %find heartbeat samples to show to the expert.
                        FindClusterExamples;

                        % Ask for expert opinion
                        [ Label bRefresh bCancel ] = ExpertUserInterface(cCentroid, cCloserExamples, cDistantExamples);

                        if(bCancel)
                            fprintf(2, 'Canceling current classification.\n');
                            if( bHaveUserInterface )                            
                                if( ishandle(UCP_struct.fig_hdl) )
                                    waitbar(0, UCP_struct.fig_hdl, 'Classification progress 0%')
                                    set(UCP_struct.Axes_hdl, 'Visible','off' );
                                end
                            end
                            break;
                        end

                        if( ~bRefresh  )

                            if( isempty(Label) )
                                bDataSkiped = true;
                                disp(['Skipping ' num2str(clust_data_m) ' heartbeats to re-cluster.'])
                            else
                                Clust_Labels(not_labeled_idx(aux_idx)) = find( strcmpi(cellstr(typical_lablist), Label) );
                                this_loop_already_labeled_idx = [ this_loop_already_labeled_idx; colvec(not_labeled_idx(aux_idx)) ];

                                hb_percent_done = (length(already_labeled_idx) + length(this_loop_already_labeled_idx))/cant_QRS_locations;
                                waitbar( hb_percent_done, UCP_struct.fig_hdl, [ 'Classification progress ' num2str(round(hb_percent_done*100)) '%' ])
                            end
                            ii = ii + 1;
                        end
                    end
                end
            end

            if( ~bCancel )

                %save the labeling done
                Labels(pending_hb_idx) = Clust_Labels;
                already_labeled_idx = [already_labeled_idx; colvec(pending_hb_idx(this_loop_already_labeled_idx))];
                pending_hb_idx = setdiff(1:cant_QRS_locations, already_labeled_idx );

                if(bHaveUserInterface)
                    if( bDataSkiped  )
                        %Maybe some fine tuning is necesary ...
                        if( ~exist('UCP_struct', 'var') || ~ishandle(UCP_struct.fig_hdl) )
                            UserControlPanel;
                        end

                        cant_hb_pending = length(pending_hb_idx);
                        fprintf(1, [num2str(cant_hb_pending) ' heartbeats remaining. Any fine tuning needed ? Change controls freely\n' ] );

                        if( cant_hb_pending < 90 ) % 10 * (k+1) * clusters
                            recommended_clusters = max(2, round(cant_hb_pending/90));
                            fprintf(2, [ 'Too few heartbeats. Recommended clusters: '  num2str(recommended_clusters) ]);
                            set(UCP_struct.CantClusters, 'Value', recommended_clusters);
                        end
                        
%                         set(UCP_struct.fig_hdl, 'Visible','on' );
                        EnableControPanel;

                        uiwait(UCP_struct.fig_hdl);
                        
                        if( ishandle(UCP_struct.fig_hdl) )
                            DisableControPanel;
%                             set(UCP_struct.fig_hdl, 'Visible','off' );
                        else
                            bCancel = true;
                            fprintf(2, 'Canceling current classification.\n');
                        end
                    end

                    %to update current user interaction.
                    ParseUserControlPanelInput;
                end
                
            end

        end

        %% Results generation

        if( bCancel )
            %Cancelation of expert input ends here.
            bCancel = false;
        else

            % Update point
            pb.checkpoint('Results generation');

            if( ~isempty(true_labels) )
                [~, this_iter_cm] = DisplayResults( 'TrueLabels', typical_lablist(true_labels,:), ...
                                                    'class_filter', class_filter, ...
                                                    'ClassifiedLabels', typical_lablist(Labels,:), ...
                                                    'ClassLabels', typical_lablist);
            else
                %not possible to build a confusion matrix without the true labels.
                this_iter_cm = nan(typical_cant_labels);
            end   
        end

        % All done
        if( bHaveUserInterface )
            if( ishandle(pb) )
                pb.checkpoint('Saving error report');
            end
        end
        
    catch MException
        %% Error handling

        if( bHaveUserInterface)
            %% with UI
            if( isempty(MException.identifier) || ~isempty(strfind(MException.identifier, 'a2hbc')) || ~isempty(strfind(MException.identifier, 'ECGwrapper'))  )
                %Our errors
                if( strfind(MException.identifier, 'ArgCheck') )
                    %Argument error, try other settings if user interface is
                    %available
                    fprintf(2, '%s', MException.message);
                    fprintf(2, '\nTry a different set of arguments with the user interface.\n');
                else
                    %other home-made errors. Make an educated exit ...
                    rethrow(MException)
                end

            else
                %% Other unknown errors -> Report it to me

                if( userNoAnswers <= maxNoAnswers && strcmpi( 'Yes', questdlg('Oh no! An unexpected error ocurred. Please take a minute to report it so I can fix it and improve A2HBC. Report error information through Internet ?', 'Unhandled error', 'Yes', 'No', 'Yes')) )
                    %report an error
                    if( ~exist('pb', 'var') )
                        % Activate the progress_struct bar.
                        pb = progress_bar(rec_filename);
                    end
                    pb.checkpoint('Saving error report');

                    timestamp = datestr(now, 'HH-MM-SS_dd-mm-yy');
                    report_filename = [tmp_path 'error_report_' timestamp '_' rec_filename '.mat' ];

                    %clean heavy variables
                    variables = whos();
                    [var_size, var_idx] = sort(cell2mat({variables(:).bytes}));
                    cum_var_size = cumsum(var_size);
                    last_var_idx = find( (cum_var_size/typical_compression_ratio) <= max_report_filesize, 1, 'last');
                    %save does not accept cell arrays as input -> workaround
                    var_names = {variables(var_idx(1:last_var_idx)).name};
                    feval( @save, report_filename, var_names{:} )

                    pb.checkpoint('Compressing');
                    
                    report_filename_gz = gzip( report_filename );
                    delete( report_filename );

                    pb.checkpoint('Reporting');

                    setpref('Internet','E_mail', 'a2hbc.errors@gmail.com');
                    setpref('Internet','SMTP_Server','smtp.gmail.com');
                    setpref('Internet','SMTP_Username','a2hbc.errors');
                    setpref('Internet','SMTP_Password', 'Any_password');
                    props = java.lang.System.getProperties;
                    props.setProperty('mail.smtp.auth','true');
                    props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
                    props.setProperty('mail.smtp.socketFactory.port','465');   

%                     local_host = 'unknown';
%                     ip_address = 'unknown';
%                     location = 'unknown';
                    
                    %just to identify reports
                    local_host = java.net.InetAddress.getLocalHost;
                    
                    try
                        ip_address = urlread('http://automation.whatismyip.com/n09230945.asp');
                    catch dummyME
                        ip_address = '';
                    end
    
                    try
                        location = urlread('http://www.ipaddresslocation.org/my-ip-address.php');
                    catch dummyME
                        location = '';
                    end

                    computer_arch = computer();
                    computer_installation = ver();
                    lcomputer_installation = length(computer_installation);
                    tab_aux = cellstr(repmat('\t\t', lcomputer_installation, 1));
                    crlf_aux = repmat('\n', lcomputer_installation, 1);
                    computer_installation = strcat( colvec({computer_installation(:).Name}), tab_aux, colvec({computer_installation(:).Version}), tab_aux, colvec({computer_installation(:).Release}) );
                    computer_installation = [char(computer_installation) crlf_aux];
                    computer_installation = cellstr(sprintf(colvec(computer_installation')));

                    sendmail({'llamedom@unizar.es' 'llamedom@electron.frba.utn.edu.ar'}, ... %recipients
                            ['[A2HBC Error reporting] ' char(local_host)], ... %subject
                            [ cellstr(location) ; cellstr(ip_address) ; cellstr(computer_arch) ; computer_installation], ... %body
                            report_filename_gz); %attach

                    pb.checkpoint('Report sent successfully');

                    delete( char(report_filename_gz) );
                    
                    fprintf(1, '\n------------------------\nReport sent successfully\n------------------------\n\nThis is the error information:\n\n');
                    
                else                
                    %try to convince user for the next time.
                    fprintf( 2, '\nOk, thank you anyway. Please consider sending a report in case this error happens again.\n\n')
                    
                    if(userNoAnswers <= maxNoAnswers)
                        userNoAnswers = userNoAnswers + 1;
                    end
                end

                rethrow(MException);
            end
        else
            %% No User interface
            
% this is obsolete now since ECGwrapper takes care about error reporting
% 
%             fprintf(2,'\n\n')
%             fprintf(2,'###########\n')
%             fprintf(2,'## ERROR ##\n')
%             fprintf(2,'###########\n')
% 
%             fprintf(2,'Recording: %s (%d/%d) \n', recording_name, this_pid, cant_pids);
% 
%             strAux = GetFunctionInvocation(mfilename, varargin);
%             
%             local_host = getenv('HOSTNAME');
%             computer_arch = computer();
% 
%             fprintf(2,'Computer: %s (%s) \n', local_host, computer_arch);
%             
%             fprintf(2, '%s\n\n', strAux);
%             
%             report = getReport(MException);
%             fprintf(2, '%s', report);
%             fprintf(2,'###########\n')
%             fprintf(2,'## ERROR ##\n')
%             fprintf(2,'###########\n')
            
            rethrow(MException)

        end
        

    end

    ConfusionMatrix(:,:,repeat_idx) = this_iter_cm;
    
    if( bHaveUserInterface && ( ~bArgCheckPassed || bInteractiveMode ) )
        %% Ask User interaction
        % only if: the arguments were incorrect or I want to use a2hbc
        % in interactive mode.
        
        UserInteraction;

    else
        %If not user interface available -> unatended mode, just exit the loop.
        %bUserExit = true;
        repeat_idx = repeat_idx + 1;
    end

end

LabelList = typical_lablist;

if( bECG_sample_provided )
    delete([ tmp_path CachedFeatMatFileName])
end

%save results when running without user interface.
if( ~bHaveUserInterface )
    
    save([tmp_path 'tmpfile_' rec_filename '_results_' op_mode '_' num2str(CantClusters) '_' num2str(iter_times) '_' num2str(cluster_presence) '_' num2str(this_pid) '.mat' ], 'Labels', 'ConfusionMatrix', 'LabelList');
    
    % flag that the program ended correctly
    setenv('A2HBC_ES', '0');
end


%%%%%%%%%%%%%%%%%%
% Other functions
%%%%%%%%%%%%%%%%%%

function posterior = CalcLDCposterior(data, ldc_struct)

[ cant_classes k ] = size(ldc_struct.mean);

posterior = nan(size(data,1), cant_classes);

for ii = 1:cant_classes
    
    centered_data = bsxfun( @minus, data, ldc_struct.mean(ii,:));

    posterior(:,ii) = exp(-0.5*sum(centered_data'.*(ldc_struct.cov(:,:,ii)*centered_data'),1)' - (ldc_struct.det(ii) + k*log(2*pi))*0.5);

end

function RunNow(obj, event_obj)

fig_hnd = gcf();
user_data = get(fig_hnd , 'User');
uiresume(fig_hnd); 

function SliderFunc(obj, event_obj)
fig_hnd = gcf();
user_data = get(fig_hnd , 'User');
NewVal = round(get(obj , 'Value'));
set(obj , 'Value', NewVal );
set(user_data.SliderLabel , 'String', [ num2str(NewVal) ' Clusters'] );

function SliderFunc2(obj, event_obj)
fig_hnd = gcf();
user_data = get(fig_hnd , 'User');
NewVal = round(get(obj , 'Value'));
set(obj , 'Value', NewVal );
strAux = [ num2str(NewVal) ' repetition'];
if(NewVal > 1)
    strAux= [strAux 's'];
end
set(user_data.SliderLabel2 , 'String', strAux );

function SliderFunc3(obj, event_obj)
fig_hnd = gcf();
user_data = get(fig_hnd , 'User');
NewVal = round(get(obj , 'Value'));
set(obj , 'Value', NewVal );
set(user_data.SliderLabel3 , 'String', [ 'Cluster majority ' num2str(NewVal) ' %'] );

function UpdateClasses(obj, event_obj)
fig_hnd = gcf();
user_data = get(fig_hnd , 'User');
Selection = get(obj , 'Value');
aux_val = get(user_data.class_filter , 'Value');
aux_str = get(user_data.class_filter , 'String');

class_filter = cellstr(aux_str(aux_val,:));

[~, aux_val] = intersect(user_data.typical_lablists{Selection}, class_filter );

set(user_data.class_filter , 'String', char(user_data.typical_lablists{Selection}) );
set(user_data.class_filter , 'Value', aux_val );


function OpModeSelect(obj, event_obj)
fig_hnd = gcf();
user_data = get(fig_hnd , 'User');
Selection = get(obj , 'Value');

if( Selection > 2 )
    set(user_data.cluster_presence , 'Enable', 'off' );
else
    set(user_data.cluster_presence , 'Enable', 'on' );
end

function FormatCheck(obj, event_obj)
fig_hnd = gcf();
user_data = get(fig_hnd , 'User');
NewVal = round(get(obj , 'Value'));
set(obj , 'Value', NewVal );
set(user_data.SliderLabel3 , 'String', [ 'Cluster majority ' num2str(NewVal) ' %'] );

function BrowseRecording(obj, event_obj)
fig_hnd = gcf();
user_data = get(fig_hnd , 'User');
recording_name = get(user_data.recording_name , 'String');
[recording_path, ~, recording_ext] = fileparts(recording_name);
[recording_name, recording_path] = uigetfile([recording_path filesep '*.*'], 'Select a recording');
if( recording_name ~= 0 )
    set(user_data.recording_name , 'String', [ recording_path recording_name]);
end

function BrowseTMPpath(obj, event_obj)
fig_hnd = gcf();
user_data = get(fig_hnd , 'User');
tmp_path = get(user_data.tmp_path , 'String');
tmp_path = uigetdir(tmp_path, 'Select temporal directory');
if( tmp_path ~= 0 )
    set(user_data.tmp_path , 'String', tmp_path);
end

function PathCheck(obj, event_obj)

tmp_path = get(obj , 'String');

if( exist( tmp_path, 'dir' ) )
    set(obj , 'String', tmp_path)
    set(obj , 'Tag', tmp_path)
else
    fprintf(2, ['Path not found. ' tmp_path '\n'] )
    set(obj , 'String', get(obj , 'Tag') )
end

function RecordingCheck(obj, event_obj)

recording = get(obj , 'String');

if( exist( recording, 'file' ) )
    set(obj , 'String', recording)
    set(obj , 'Tag', recording)
else
    fprintf(2, ['File not found. ' recording '\n'] )
    set(obj , 'String', get(obj , 'Tag') )
end

function CheckNumeric(obj, event_obj)

number = get(obj , 'String');

if( isnan(str2double(number)) )
    fprintf(2, ['Not a valid numeric value: ' number '\n'] )
    set(obj , 'String', get(obj , 'Tag') )
else
    set(obj , 'String', number)
    set(obj , 'Tag', number)
end

