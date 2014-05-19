function [payload, global_struct ] = ECG_avg( ECG, ECG_start_end_idx, ECG_header, ECG_annotations, global_struct)

ECG_leads = size(ECG, 2);
ECG_length = diff(ECG_start_end_idx) + 1;

% do some checksum
if( isempty(ECG_annotations) )
    payload = [diff(ECG_start_end_idx)+1, rem(sum(ECG(ECG_start_end_idx(1):ECG_start_end_idx(2),:)),1024) ];
else
    payload = [size(ECG_annotations.time,1), rem(sum(ECG(ECG_annotations.time,:)),1024) ];
end

if( isfield(global_struct, 'ECGmean' ) )
    aux_sum = global_struct.ECGsize + ECG_length;
    global_struct.ECGmean = 1/aux_sum*(global_struct.ECGsize .* global_struct.ECGmean + ECG_length .* mean(ECG(ECG_start_end_idx(1):ECG_start_end_idx(2),:)));
    global_struct.ECGsize = aux_sum;
else
    global_struct.ECGmean = mean(ECG(ECG_start_end_idx(1):ECG_start_end_idx(2),:));
    global_struct.ECGsize = ECG_length;
end

% fprintf(1, ['ECG mean: [' repmat( '%5.4e ', 1, ECG_leads) ']\n'], global_struct.ECGmean);

