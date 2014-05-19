function qrs_sample_idx = plot_n_QRS(ECG, t_idx, subs_win_size_samples, FP_heartbeats_hdl, cant_qrs )

cla(FP_heartbeats_hdl);
hold(FP_heartbeats_hdl, 'on')
lt_idx = length(t_idx);
cant_qrs = min(lt_idx, cant_qrs);
t_idx = rowvec(t_idx);

qrs_sample_idx = randsample(lt_idx, cant_qrs );

hb_idx = arrayfun(@(a,b)(colvec(a:b)), t_idx(qrs_sample_idx) - 4*subs_win_size_samples, ...
                                       t_idx(qrs_sample_idx) + 4*subs_win_size_samples, ...
                                       'UniformOutput', false);
hb_idx = cell2mat(hb_idx);

ECG = reshape(ECG(hb_idx,1), size(hb_idx));
ecg_offset = max(max(ECG)) - min(min(ECG));
plot(FP_heartbeats_hdl, ECG - repmat((0:cant_qrs-1)*ecg_offset,size(hb_idx,1),1) );
plot(FP_heartbeats_hdl, 4*[subs_win_size_samples subs_win_size_samples], [ 1 -(cant_qrs-1)]*ecg_offset, 'r--' );

axes(FP_heartbeats_hdl);
for ii = 1:length(qrs_sample_idx)
    text( -20, -(ii-1)*ecg_offset, num2str(t_idx(qrs_sample_idx(ii))) );    
end

xlim(FP_heartbeats_hdl, [-30 size(ECG,1)+20]);
hold(FP_heartbeats_hdl, 'off');
