
%centroid idx is the first
clust_sample_centroid_idx = clust_dist2mu_sorted_idx(1);

centroid_FM_idx = not_labeled_idx(aux_idx(aux_idx2(clust_sample_centroid_idx)));

centroid_QRS_idx = (centroid_FM_idx);

aux_hb_seq =    max(1, QRS_locations(centroid_QRS_idx) - round(2*sampling_rate)): ...
                min(ECG_header.nsamp, QRS_locations(centroid_QRS_idx) + round(2*sampling_rate));

aux_hb_seq2 =   max(1, centroid_FM_idx - 5):min(cant_QRS_locations, centroid_FM_idx + 2);

aux_hb_seq3 =    max(1, QRS_locations(centroid_QRS_idx) - round(0.2*sampling_rate)): ...
                min(ECG_header.nsamp, QRS_locations(centroid_QRS_idx) + round(0.5*sampling_rate));


cCentroid = {   (QRS_locations(centroid_QRS_idx) - aux_hb_seq(1) + 1) ... %relative position of the heartbeat
                ECG_total(aux_hb_seq,:) ... % ecg excerpt
                centroid_FM_idx - aux_hb_seq2(1) + 1 ... % relative position of the heartbeat 
                RR_intervals(aux_hb_seq2)/sampling_rate ... % RR interval sequence
                (QRS_locations(centroid_QRS_idx) - aux_hb_seq3(1) + 1) ... %relative position of the heartbeat                
                ECG_total(aux_hb_seq3,:) ... % morph excerpt
                clust_data_m ... % cluster size
                };
            
%closer samples
quartile_size = round(clust_data_m * 0.25);
if( quartile_size > CantSamplesClose )
    clust_samples_close_idx = clust_dist2mu_sorted_idx( 1 + randsample(quartile_size, CantSamplesClose ) );
else
    if( clust_data_m < 2 )
        % the centroid itself
        clust_samples_close_idx = clust_dist2mu_sorted_idx;
    else
        clust_samples_close_idx = clust_dist2mu_sorted_idx( 2:(2+quartile_size-1) );
    end
end

clust_samples_close_idx = sort(clust_samples_close_idx);
clust_samples_close_FM_idx = not_labeled_idx(aux_idx(aux_idx2(clust_samples_close_idx)));
clust_samples_close_QRS_idx = (clust_samples_close_FM_idx);


aux_ECGslices = [];
aux_RRslices = [];
aux_MorphSlices = [];
aux_rel_posQRS = [];
aux_rel_posRR = [];
aux_rel_MorphposQRS = [];

   
    aux_seq = colvec(1:length(clust_samples_close_QRS_idx));
    
    aux_hb_seq = arrayfun(@(a)( max(1, QRS_locations(clust_samples_close_QRS_idx((a))) - round(2*sampling_rate)): ...
                                min(ECG_header.nsamp, QRS_locations(clust_samples_close_QRS_idx((a))) + round(2*sampling_rate)) ) , ...
                           aux_seq, 'UniformOutput', false);

    aux_hb_seq2 = arrayfun(@(a)(    max(1, clust_samples_close_FM_idx((a)) - 5): ...
                                    min(cant_QRS_locations, clust_samples_close_FM_idx((a)) + 2) ) , ...
                           aux_seq, 'UniformOutput', false);

    aux_hb_seq3 = arrayfun(@(a)( max(1, QRS_locations(clust_samples_close_QRS_idx((a))) - round(0.2*sampling_rate)): ...
                                min(ECG_header.nsamp, QRS_locations(clust_samples_close_QRS_idx((a))) + round(0.5*sampling_rate)) ) , ...
                           aux_seq, 'UniformOutput', false);

    %ECG slices
    aux_ECGslices = [aux_ECGslices; cellfun(@(a)(ECG_total(a,:)), aux_hb_seq, 'UniformOutput', false)];
    %RR slices
    aux_RRslices = [aux_RRslices; cellfun(@(a)(RR_intervals(a)/sampling_rate), aux_hb_seq2, 'UniformOutput', false)];
    %Morph slices
    aux_MorphSlices = [aux_MorphSlices; cellfun(@(a)(ECG_total(a,:)), aux_hb_seq3, 'UniformOutput', false)];
    
    %relative positions
    aux_rel_posQRS = [aux_rel_posQRS; cellfun(@(a,b)(QRS_locations(clust_samples_close_QRS_idx((a))) - b(1) + 1), num2cell(aux_seq), aux_hb_seq)];
    aux_rel_posRR =  [aux_rel_posRR; cellfun(@(a,b)(clust_samples_close_FM_idx((a)) - b(1) + 1), num2cell(aux_seq), aux_hb_seq2)];
    aux_rel_MorphposQRS = [aux_rel_MorphposQRS; cellfun(@(a,b)(QRS_locations(clust_samples_close_QRS_idx((a))) - b(1) + 1), num2cell(aux_seq), aux_hb_seq3)];
                       

if( iscell(aux_ECGslices) )
    cCloserExamples = { aux_rel_posQRS ... %relative position of the heartbeat
                        aux_ECGslices ... % ecg excerpt
                        aux_rel_posRR ... % relative position of the heartbeat 
                        aux_RRslices ... % RR interval sequence
                        aux_rel_MorphposQRS ... % Morph relative position of the heartbeat 
                        aux_MorphSlices ... % Morph interval sequence
                        };
else
    cCloserExamples = { aux_rel_posQRS ... %relative position of the heartbeat
                        {aux_ECGslices} ... % ecg excerpt
                        aux_rel_posRR ... % relative position of the heartbeat 
                        {aux_RRslices} ... % RR interval sequence
                        aux_rel_MorphposQRS ... % Morph relative position of the heartbeat 
                        {aux_MorphSlices} ... % Morph interval sequence
                        };
    
end

%distant samples
if( quartile_size > CantSamplesFar )
    clust_samples_far_idx = clust_dist2mu_sorted_idx( clust_data_m - randsample(quartile_size, CantSamplesFar ) + 1  );
else
    if( clust_data_m < 3 )
        % all examples
        clust_samples_far_idx = clust_dist2mu_sorted_idx;
    else
        clust_samples_far_idx = clust_dist2mu_sorted_idx(3*quartile_size:clust_data_m);
    end
end


clust_samples_far_idx = sort(clust_samples_far_idx);
clust_samples_far_FM_idx = not_labeled_idx(aux_idx(aux_idx2(clust_samples_far_idx)));
clust_samples_far_QRS_idx = (clust_samples_far_FM_idx);

aux_ECGslices = [];
aux_RRslices = [];
aux_MorphSlices = [];
aux_rel_posQRS = [];
aux_rel_posRR = [];
aux_rel_MorphposQRS = [];

    
    aux_seq = colvec(1:length(clust_samples_far_QRS_idx));
    
    aux_hb_seq = arrayfun(@(a)( max(1, QRS_locations(clust_samples_far_QRS_idx((a))) - round(2*sampling_rate)): ...
                                min(ECG_header.nsamp, QRS_locations(clust_samples_far_QRS_idx((a))) + round(2*sampling_rate)) ) , ...
                           aux_seq, 'UniformOutput', false);

    aux_hb_seq2 = arrayfun(@(a)(    max(1, clust_samples_far_FM_idx((a)) - 5): ...
                                    min(cant_QRS_locations, clust_samples_far_FM_idx((a)) + 2) ) , ...
                           aux_seq, 'UniformOutput', false);

    aux_hb_seq3 = arrayfun(@(a)( max(1, QRS_locations(clust_samples_far_QRS_idx((a))) - round(0.2*sampling_rate)): ...
                                min(ECG_header.nsamp, QRS_locations(clust_samples_far_QRS_idx((a))) + round(0.5*sampling_rate)) ) , ...
                           aux_seq, 'UniformOutput', false);
                       
    %ECG slices
    aux_ECGslices = [aux_ECGslices; cellfun(@(a)(ECG_total(a,:)), aux_hb_seq, 'UniformOutput', false)];
    %RR slices
    aux_RRslices = [aux_RRslices; cellfun(@(a)(RR_intervals(a)/sampling_rate), aux_hb_seq2, 'UniformOutput', false)];
    %Morph slices
    aux_MorphSlices = [aux_MorphSlices; cellfun(@(a)(ECG_total(a,:)), aux_hb_seq3, 'UniformOutput', false)];
    
    %relative positions
    aux_rel_posQRS = [aux_rel_posQRS; cellfun(@(a,b)(QRS_locations(clust_samples_far_QRS_idx((a))) - b(1) + 1), num2cell(aux_seq), aux_hb_seq)];
    aux_rel_posRR =  [aux_rel_posRR; cellfun(@(a,b)(clust_samples_far_FM_idx((a)) - b(1) + 1), num2cell(aux_seq), aux_hb_seq2)];
    aux_rel_MorphposQRS = [aux_rel_MorphposQRS; cellfun(@(a,b)(QRS_locations(clust_samples_far_QRS_idx((a))) - b(1) + 1), num2cell(aux_seq), aux_hb_seq3)];


if( iscell(aux_ECGslices) )
    cDistantExamples ={ aux_rel_posQRS ... %relative position of the heartbeat
                        aux_ECGslices ... % ecg excerpt
                        aux_rel_posRR ... % relative position of the heartbeat 
                        aux_RRslices ... % RR interval sequence
                        aux_rel_MorphposQRS ... % Morph relative position of the heartbeat 
                        aux_MorphSlices ... % Morph interval sequence
                        };
else
    %force cell contents.
    cDistantExamples ={ aux_rel_posQRS ... %relative position of the heartbeat
                        {aux_ECGslices} ... % ecg excerpt
                        aux_rel_posRR ... % relative position of the heartbeat 
                        {aux_RRslices} ... % RR interval sequence
                        aux_rel_MorphposQRS ... % Morph relative position of the heartbeat 
                        {aux_MorphSlices} ... % Morph interval sequence
                        };
end

clear aux

