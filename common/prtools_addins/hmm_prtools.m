%HMM Hidden Markov Model adapted for PRtools
%
%   [W,R,S,M] = QDC(A,R,S,M)
%   W = A*QDC([],R,S)
%
% INPUT
%   A    Dataset
%   R,S	 Regularization parameters, 0 <= R,S <= 1 
%        (optional; default: no regularization, i.e. R,S = 0)
%   M    Dimension of subspace structure in covariance matrix (default: K,
%        all dimensions)
%
% OUTPUT
%   W    Quadratic Bayes Normal Classifier mapping
%   R    Value of regularization parameter R as used 
%   S    Value of regularization parameter S as used
%   M    Value of regularization parameter M as used
%
% DESCRIPTION
% Computation of the quadratic classifier between the classes of the dataset
% A assuming normal densities. R and S (0 <= R,S <= 1) are regularization
% parameters used for finding the covariance matrix by
% 
%   G = (1-R-S)*G + R*diag(diag(G)) + S*mean(diag(G))*eye(size(G,1))
%
% This covariance matrix is then decomposed as G = W*W' + sigma^2 * eye(K),
% where W is a K x M matrix containing the M leading principal components
% and sigma^2 is the mean of the K-M smallest eigenvalues.
%
% 
% 
% The use of soft labels is supported. The classification A*W is computed by
% NORMAL_MAP.
%
% If R, S or M is NaN the regularisation parameter is optimised by REGOPTC.
% The best result are usually obtained by R = 0, S = NaN, M = [], or by
% R = 0, S = 0, M = NaN (which is for problems of moderate or low dimensionality
% faster). If no regularisation is supplied a pseudo-inverse of the
% covariance matrix is used in case it is close to singular.
%
% EXAMPLES
% See PREX_MCPLOT, PREX_PLOTC.
%
% REFERENCES
% 1. R.O. Duda, P.E. Hart, and D.G. Stork, Pattern classification, 2nd
% edition, John Wiley and Sons, New York, 2001. 
% 2. A. Webb, Statistical Pattern Recognition, John Wiley & Sons, 
% New York, 2002.
%
% SEE ALSO
% MAPPINGS, DATASETS, REGOPTC, NMC, NMSC, LDC, UDC, QUADRC, NORMAL_MAP

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: qdc.m,v 1.7 2008/03/20 09:25:10 duin Exp $

function result_out = hmm_prtools(varargin)

	prtrace(mfilename);
    
    cEmissionCalcMode = { 'normal' 'perPreviousClass' 'custom'};
    
%     'normal': One classifier per state or hidden node (e.g. one cov marix and mean vector per state)
%     'perPreviousClass': For a C class classification problem, one classifier given that class i happened (e.g. Given one Supraventricular beat, the same cov matrix for all C classes, and one mean vector for each class)
%     'custom': Custom calculation of cov matrices and mean vectors for each state.
    
    %input parsing
    p = inputParser;   % Create instance of inputParser class.

    p.addOptional('dsAny',    [], @(x)(isdataset(x) ) );
    p.addOptional('wMapping', [], @(x)(ismapping(x) ) );
    p.addParamValue('CantPrevBeats', 0, @(x)(isnumeric(x) && x >= 0 ) );
    p.addParamValue('TransMat', [], @(x)(isnumeric(x) && all(all(x >= 0)) ) );
    p.addParamValue('InitialStateProb', [], @(x)(isnumeric(x) && all(x >= 0) ) );
    p.addParamValue('EmmissionMapping', weighted_ldc([],true), @(x)( ismapping(x) ) );
    p.addParamValue('TiedStates', [], @(x)(isnumeric(x)) );
    p.addParamValue('Details', false, @(x)(islogical(x)) );
    p.addParamValue('PooledCovMat', false, @(x)(islogical(x)) );
    p.addParamValue('EmissionEstimationMode', 'normal', @(x)any(strcmpi(x,cEmissionCalcMode)) );

    p.parse( varargin{:} );

    dsAny = p.Results.dsAny;
    wMapping = p.Results.wMapping;
    CantPrevBeats = p.Results.CantPrevBeats;
    TransMat = p.Results.TransMat;
    InitialStateProb = p.Results.InitialStateProb;
    wEmmissionMapping = p.Results.EmmissionMapping;
    TiedStates = p.Results.TiedStates;
    EmissionEstimationMode = p.Results.EmissionEstimationMode;
    bDetails = p.Results.Details;
    bPooledCovMat = p.Results.PooledCovMat;
    
	if( nargin < 1 || isempty(dsAny) )       
        %% No input arguments: 
		
        result_out = mapping(mfilename,varargin); % return an untrained mapping.
        result_out = setname(result_out,'HMM');
	
    elseif( ~isempty(wMapping) && istrained(wMapping) && isdataset(dsAny) ) 
        %% mapping
        
		m = getsize(dsAny,1);
		[k c] = getsize(wMapping);
        
        result_out = zeros(m,c);
        
        [ obs_labels, ~, dsAny ] = ds_convert_prtools_pmtk(dsAny);
        
        recording_index = getident(dsAny, 'recording_index');

        % hare una observacion de cada registro
        [recordings, ~, rec_idx] = unique(recording_index);

        cant_recordings = length(recordings);

        model = +wMapping;

        predicted_seq = cell(cant_recordings,1);
        lStateGroups = length(model.TiedStatesGroups);

        for ii = 1:cant_recordings

            each_rec_idx = find( ii == rec_idx );
            
            logB = nan(model.hmm_states, length(each_rec_idx) );
            
            %protect from NaNs
            aux = +dsAny(each_rec_idx,:);
            aux_mean = nanmean(aux);
            aux_nan_idx = find(isnan(aux));
            [aux_nan_row aux_nan_col] = ind2sub( size(aux), aux_nan_idx);
            
            for jj = 1:length(aux_nan_row)
                aux(aux_nan_row(jj), aux_nan_col(jj)) = aux_mean(aux_nan_col(jj));
            end
            
            dsAny(each_rec_idx,:) = aux;
            clear aux
            %END protect from NaNs
            
            for jj = 1:lStateGroups

                this_states = find(model.TiedStatesGroups_idx == jj);
                
                switch(model.EmissionEstimationMode)

                    case 'normal' 

                        logB = dsAny(each_rec_idx,:) * model.emission;
%                         logB = normalize( (+logB)', 1);
                        logB = (+logB)';

                    case 'perPreviousClass' 

                        aux = dsAny(each_rec_idx,:) * model.emission{jj};
%                         logB(this_states,:) = normalize( (+aux)', 1);
                        logB(this_states,:) = (+aux)';

                    case 'custom' 

                        aux = dsAny(each_rec_idx,:) * model.emission{jj};
%                         aux = normalize( (+aux)', 1);
                        aux = (+aux)';
                        tied_lobB = size(aux, 1);
                        tied_classes = length(this_states);
                        
                        if( tied_lobB == tied_classes )
                            logB(this_states,:) = aux; 
                        else
                            logB(this_states,:) = repmat(aux, tied_classes/tied_lobB, 1);
                        end

                    otherwise 

                        error('Unrecognized case.');                
                end

            end
            
            [predicted_seq{ii}, ~, ~] = hmmViterbiC(log(model.pi+eps), log(model.A+eps), log(logB));
%             [predicted_seq{ii}, ~, ~] = hmmViterbiCM(log(model.pi+eps), log(model.A+eps), log(logB));
            
% %             conv_idx = [4 1 2 4 3];
%             figure(99)
% %             plot(conv_idx(obs_labels{ii}), 'bo-' )
%             plot(obs_labels{ii}, 'bo-' )
%             hold on
%             ps = unchange_labels(predicted_seq{ii}, model.CantPrevBeats+1, c);
%             ps2 = unchange_labels(predicted_seq2{ii}, model.CantPrevBeats+1, c);
%             plot(ps, 'rx-' )
%             plot(ps2, 'gs-' )
%             hold off
            
        end
        
        predicted_seq = cell2mat(predicted_seq')';
        predicted_changed = predicted_seq;
        predicted_seq = unchange_labels(predicted_changed, model.CantPrevBeats+1, c);
        
        if(model.bDetails)
            
            true_labels = getnlab(dsAny);
            
            %take class strings initials letters
            learned_lablist = getlabels(wMapping);
            
            test_lablist = getlablist(dsAny);
            
            orig_lablist = (learned_lablist(:,1));
            
            lablist = orig_lablist';
            beat_sequences = orig_lablist;
            cant_labels = length(lablist);

            [~, test_lablist_idx, orig_lablist_idx] = intersect(cellstr(test_lablist), cellstr(learned_lablist) );

            test_cant_labels = size(test_lablist,1);
            
            conversion_idx = repmat(c+1,test_cant_labels,1);
            
            if( ~isempty(test_lablist_idx) )
                conversion_idx(test_lablist_idx) = orig_lablist_idx;
            end
            
            if( test_cant_labels > c )
                conv_lablist = [ orig_lablist; '~'];
%                 warning('More classes than trained in evaluation dataset, "~" means not learned classes.')
            else
                conv_lablist = orig_lablist;
            end
            conv_lablist = conv_lablist(:);
            
            % descarto las clases no aprendidas
            true_labels = conversion_idx(true_labels);
            learned_examples_idx = find(true_labels<=c);
            
            %Generate those labels for each state
            dsAux = setlabels(dsAny, true_labels);
            [dsAux selected_idx] = seldat(dsAux, 1:3);
            test_true_labels = ds_convert_prtools_pmtk(dsAux);
            true_states_labels = change_labels(test_true_labels, model.CantPrevBeats+1, c);
            true_states_labels = cell2mat(true_states_labels')';
            
            true_labels = true_labels(learned_examples_idx);
            predicted_seq_str = conv_lablist(predicted_seq(learned_examples_idx));
            true_labels_str = conv_lablist( true_labels );
            
            %Generate a lablist of strings for each hidden state
            for ii = 2:(model.CantPrevBeats+1)
                aux = repmat(lablist,size(beat_sequences,1),1);
                beat_sequences = [ aux(:)  repmat(beat_sequences,cant_labels,1) ];
            end
            lablist = beat_sequences;
            
            true_states_labels_str = lablist(true_states_labels,:);
            predicted_changed = lablist(predicted_changed(selected_idx),:);

            %Displays confmat of hidden states
            confmat(true_states_labels_str, predicted_changed);

            %Generate filters to group by previous class sequence history,
            %and to measure the relative importance of each group in the
            %evaluation dataset.
            cant_labels_in = [];
            labels_in = true(size(true_labels));
            jj = 1;

            for ii = 1:c:model.hmm_states
                labels_in(:,jj) = true_states_labels >= ii & true_states_labels < (ii+c);
                cant_labels_in(jj) = sum(labels_in(:,jj));
                jj = jj + 1;
            end

            [~, importance_idx] = sort(cant_labels_in, 'descend');

            c_aux = [];
            jj = 1;
            for ii = importance_idx
                prev_beats(jj,:) = true_states_labels_str( find(labels_in(:,ii), 1,'first') , 1:model.CantPrevBeats );
                disp(prev_beats(jj,:))
                confmat(true_labels_str(labels_in(:,ii)), predicted_seq_str(labels_in(:,ii)) );
                [ConfusionMat,NE,lab_list_en_register_out,lab_list_en_trained_classifier ] = confmat(true_labels_str(labels_in(:,ii)), predicted_seq_str(labels_in(:,ii)) );

                % Tengo que compatibilizar la matriz de confusion resultante ya que
                % habra distintas clases en register_out y trained_classifier.
                [NLAB11,NLAB12] = renumlab(orig_lablist, lab_list_en_register_out );
                [NLAB21,NLAB22] = renumlab(orig_lablist, lab_list_en_trained_classifier );

                cc_aux = zeros(c);

                cc_aux(NLAB12,NLAB22) = ConfusionMat;

                c_aux(:,:,jj) = cc_aux;

                jj = jj + 1;
            end

            c_aux = [ c_aux; sum(c_aux) ];
            c_aux = [ c_aux sum(c_aux,2) ];

            for ii = 1:c
                disp(orig_lablist(ii))
                importance_i = squeeze( round( c_aux(ii,end,:) * 1/sum(c_aux(ii,end,:)) * 100));
                sens_i = squeeze(round( c_aux(ii,ii,:) ./ c_aux(ii,end,:) * 100));
                pospred_i = squeeze(round(c_aux(ii,ii,:) ./ c_aux(end,ii,:) * 100));

                aux = [importance_i sens_i pospred_i]';
                [dummy, importance_idx] = sort(importance_i', 'descend');
                aux_str = [prev_beats(importance_idx,:) repmat('   ', size(prev_beats,1) , 1)]';
                fprintf(1, '           %s\n', aux_str(:) )
                disp([['importance '; ...
                       'sensitiv   '; ...
                       'pospred    '] num2str(aux(:,importance_idx)) ] );

                fprintf(1, '\n' )
                fprintf(1, 'Total Se: %d\n', round( sum(c_aux(ii,ii,:)) / sum(c_aux(ii,end,:)) * 100) )
                fprintf(1, 'Total +P: %d\n\n', round( sum(c_aux(ii,ii,:)) / sum(c_aux(end,ii,:)) * 100) )

            end
            
            
        end
                
        %Simulates a soft output, given that viterbi algorithm does not
        %output probabilities.
        result_out( sub2ind(size(result_out), 1:size(result_out,1), predicted_seq' ) ) = 1;
        
        %build the result dataset.
        result_out = setdata(dsAny, result_out, getlabels(wMapping));
        
        
    elseif( isdataset(dsAny) ) 
        %% training
		
		islabtype(dsAny,'crisp'); % Assert A has the right labtype.
		isvaldfile(dsAny,2,2); % at least 2 objects per class, 2 classes

		[m,k,c] = getsize(dsAny);
        
        hmm_states = c ^ (CantPrevBeats+1);
        
        model.hmm_states = hmm_states;
        model.CantPrevBeats = CantPrevBeats;
        model.bDetails = bDetails;
        
        [obs_labels_train, dummy, dsAny] = ds_convert_prtools_pmtk(dsAny);
        
        obs_labels_train = change_labels(obs_labels_train, CantPrevBeats+1, c);
        
        obs_labels_train_stacked = cell2mat(obs_labels_train')';
        
        %% transition matrix
        if( isempty(TransMat) )
            TransMat = countTransitions(obs_labels_train, hmm_states );
        end
        
        if( bDetails )
            disp('Estimated Transition Matrix')
            disp(TransMat)
        end
        
        TransMat(TransMat == 0) = 1;
%         TransMat(:,1:c:hmm_states) = 0.5*TransMat(:,1:c:hmm_states);
%         TransMat(:,2:c:hmm_states) = 1.5*TransMat(:,2:c:hmm_states);
%         TransMat(:,3:c:hmm_states) = 2*TransMat(:,3:c:hmm_states);

        model.A  = normalize(TransMat, 2);
        
        %% initial state prior
        if( isempty(InitialStateProb) )
            seqidx   = cumsum([1, cellfun(@(seq)size(seq, 2), obs_labels_train')]);
            InitialStateProb = rowvec(histc(obs_labels_train_stacked(seqidx(1:end-1)), 1:hmm_states));
        end
        
        if( bDetails )
            disp('Estimated Initial State Probability')
            disp(InitialStateProb)
        end
        
        model.pi = normalize(InitialStateProb);
        
        %% state emission distributions estimation
        
        switch(EmissionEstimationMode)
            
            case 'normal' 
                
                TiedStates = ones(hmm_states,1);
                [StateGroups, dummy, StateGroups_idx ] = unique(TiedStates);
                train_lablist = getlablist(dsAny);
                
                aux_prior = getprior(dsAny);
                new_prior = nan(1,hmm_states);
                for ii = 1:c
                    new_prior(ii:c:hmm_states) = aux_prior(ii);
                end
                
                %Add new labels, one per state, to estimate one emission
                %per state.
                dsAny = addlabels(dsAny, obs_labels_train_stacked, 'dummy_lab');        
                
                new_prior = new_prior *1/(sum(new_prior));
                dsAny = setprior(dsAny, new_prior);

            case 'perPreviousClass' 
                
                TiedStates = repmat(1:(c^CantPrevBeats),c,1);
                TiedStates = (TiedStates(:))';
                [StateGroups, ~, StateGroups_idx ] = unique(TiedStates);
                
            case 'custom' 
                
                %Check for TiedStates
                if( ~isempty(TiedStates) )
                    if( length(TiedStates) ~= hmm_states )
                        error(['TiedStates is a numeric vector of length ' num2str(hmm_states) ])
                    end

                    [StateGroups, ~, StateGroups_idx ] = unique(TiedStates);

                else
                    warning( 'TiedStates option not defined. Assuming "full" mode.' )
                    StateGroups = 1:hmm_states;
                    StateGroups_idx = StateGroups;
                end
                
            otherwise 
                error('Unrecognized case.');
        end
        
        model.TiedStatesGroups = StateGroups;
        model.TiedStatesGroups_idx = StateGroups_idx;
        
        model.EmissionEstimationMode = EmissionEstimationMode;
        
        lStateGroups = length(StateGroups);
 
        model.emission = cell(lStateGroups,1);

        if( bDetails )
            disp('Estimated Emisssion Probabilities')
            disp('---------------------------------')
        end
        
        if( bPooledCovMat )
            wAux = weighted_ldc( dsAny, true);
            wAux = +wAux;
            PooledCovMat = wAux.cov;
            PooledCovMatDets = wAux.det;
        end
        
        for ii = 1:lStateGroups

            if( strcmpi(EmissionEstimationMode, 'normal') )
                examples_per_state_filter = true(size(obs_labels_train_stacked));
            else
                examples_per_state_filter = false(size(obs_labels_train_stacked));

                this_states = rowvec(find(StateGroups_idx == ii));

                for jj = this_states
                    examples_per_state_filter = examples_per_state_filter | obs_labels_train_stacked == jj;
                end
                
            end
            
            %Actual state emission estimation
            wMapping = dsAny(examples_per_state_filter,:) * wEmmissionMapping;
            
            if( bPooledCovMat )
                if( strcmpi(EmissionEstimationMode, 'custom') )
                    emis_aux = +wMapping;
                    emis_aux.cov = PooledCovMat;
                    emis_aux.det = PooledCovMatDets;
                    wMapping = setdata(wMapping, emis_aux);
                else
                    warning('Ignoring PooledCovMat parameter');
                end
            end
            
            if( bDetails )
                disp(['State ' num2str(ii) ])
                disp('Data: ' )
                dsAny(examples_per_state_filter,:)
                emis_aux = +wMapping;
                disp('Emission mean: ' )
                disp(emis_aux.mean)
                disp('Emission cov: ' )
                disp(emis_aux.cov)
            end
            
            switch(EmissionEstimationMode)
            
                case 'normal' 

                    model.emission = wMapping;

                case 'perPreviousClass' 
                    
                    model.emission{ii} = wMapping;

                case 'custom' 

                    model.emission{ii} = wMapping;

                otherwise 
                    
                    error('Unrecognized case.');                
            end
            
        end
        
        if( strcmpi(EmissionEstimationMode, 'normal') ) 
            result_out = mapping(mfilename, 'trained', model, train_lablist, k, c);
        else
            result_out = mapping(mfilename, 'trained', model, getlablist(dsAny), k, c);
        end
        
    else
        error('Illegal call')     % this should not happen
	end


return;



function obs_labels_train = change_labels(obs_labels_train, hmm_mem, cant_clases)

%asumo el mismo label a los hmm_mem primeros latidos.
cant_patients = length(obs_labels_train);
for ii = 1:cant_patients
    aux = obs_labels_train{ii};
    obs_labels_train{ii} = [repmat(aux(1),1,hmm_mem) aux];
end

first_idx = cumsum([1; cellfun( @length, obs_labels_train)]);
last_idx = cumsum([0; cellfun( @length, obs_labels_train)]);
first_idx = first_idx(1:end-1);
last_idx = last_idx(2:end);


obs_labels_train_orig = cell2mat(obs_labels_train')';
obs_labels_train = nan(size(obs_labels_train_orig));


all_labels = (1:cant_clases)';
cant_labels = length(all_labels);
beat_sequences = all_labels;

for ii = 2:hmm_mem
%     beat_sequences = [ [repmat(1,3,1);repmat(2,3,1);repmat(3,3,1)] repmat((1:3)',3,1) ];
    aux = repmat(all_labels',size(beat_sequences,1),1);
    beat_sequences = [ aux(:)  repmat(beat_sequences,cant_labels,1) ];
end
   
seq_length = size(beat_sequences,2);
first_beats = repmat(first_idx(:)',hmm_mem,1) + repmat((0:(hmm_mem-1))', 1, size(first_idx,1));
first_beats = first_beats(:);
ii = 1;

for jj = 1:size(beat_sequences,1)
    beat_sequence = beat_sequences(jj,:);
    Beats_idx = find(obs_labels_train_orig == beat_sequence(seq_length));
    Beats_idx = setdiff(Beats_idx, first_beats);
    bAux = true(length(Beats_idx),1);
    
    for kk = 1:seq_length-1
        bAux = bAux & obs_labels_train_orig(Beats_idx-kk) == beat_sequence(seq_length-kk);
    end    
    Beats_idx = Beats_idx(bAux);
    obs_labels_train(Beats_idx) = ii;
    ii = ii + 1;
end

cant_patients = length(first_idx);
aux = cell(cant_patients,1);
for ii = 1:cant_patients
    aux{ii} = obs_labels_train(first_idx(ii)+hmm_mem: last_idx(ii))';
end
obs_labels_train = aux;


function new_lab = unchange_labels(old_lab, hmm_mem, cant_clases)

new_lab = nan(size(old_lab));

for ii = 1:cant_clases
    for jj = ii:cant_clases:cant_clases^hmm_mem
        new_lab(jj == old_lab) = ii;
    end
end



function [ obs_labels, observations, dsAny ] = ds_convert_prtools_pmtk(dsAny)

labels = getnlab(dsAny);
data = +dsAny;

% patient_index = getident(dsAny, 'patient_index');
recording_index = getident(dsAny, 'recording_index');
beat_time = getident(dsAny, 'beat_time');

user = getuser(dsAny);
ResultsFeatlab = getfeatlab(dsAny);


% hare una observacion de cada registro
[recordings, dummy, rec_idx] = unique(recording_index);

cant_recordings = length(recordings);

observations = cell(cant_recordings,1);
obs_labels = cell(cant_recordings,1);

for ii = 1:cant_recordings
    
    each_rec_idx = find( ii == rec_idx );

%     bNanObservations = any(isnan(data(each_rec_idx,:)),2);
%     
%     %para no perder la sequencia cambiamos los nan por las medianas en cada feature por cada registro.
%     if( sum(bNanObservations) > 0 )
%        
%         this_median = nanmedian(data(each_rec_idx,:));
%        
%         for jj = find(any(isnan(data(each_rec_idx,:))))
%             thisNanObservations_idx = find(isnan(data(each_rec_idx,jj)));
%             data(each_rec_idx(thisNanObservations_idx),jj) = this_median(jj);
%         end
%         
%     end
    
    %la ordeno temporalmente
    [dummy, time_idx] = sort(beat_time(each_rec_idx));
    
    observations{ii} = data(each_rec_idx(time_idx),:)';
    obs_labels{ii} = labels(each_rec_idx(time_idx),:)';
   
    if( nargout > 2 )
        %saco el dataset ordenado temporalmente.
        dsAny(each_rec_idx,:) = dsAny(each_rec_idx(time_idx),:);
    end
    
end
