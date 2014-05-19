function [realizations, aux_idx, seq_idx] = pack_signal(signal, anns_refined, realization_limit, substract_mean)
% Description:
% Packs several realizations synchronized at each sample in anns_refined
% (QRS detections), in a window from anns_refined(i) - realization_limit(1)
% to anns_refined(i) + realization_limit(2)
% 
% Returns: realizations is a signal pack where dimensions obtained with
% size indicates:
% [ ensemble_size nsig cant_anns ] = size(realizations);
% 

if( nargin < 4)
    substract_mean = false;
end

if( isnumeric(realization_limit) && length(realization_limit) == 1 )
    realization_limit = round([realization_limit/2 realization_limit/2]);
end

[nsamp, nsig] = size(signal);

cant_anns = length(anns_refined);
seq_idx = find( (anns_refined - realization_limit(1)) >= 0 & (anns_refined + realization_limit(2)) <= nsamp );

if( length(seq_idx) ~= cant_anns )
    warning('Some annotations discarded since they were beyond the signal limits.')
    cant_anns = length(seq_idx);
end

ensemble_size = sum(realization_limit);


aux_idx = arrayfun(@(a)((anns_refined(a) - realization_limit(1) + 1):anns_refined(a) + realization_limit(2)), ...
                        seq_idx, 'UniformOutput', false);

if(substract_mean)                    
    realizations = cellfun(@(a)( bsxfun( @minus, signal(a,:), mean(signal(a,:))) ), aux_idx, 'UniformOutput', false);
else
    realizations = cellfun(@(a)( signal(a,:) ), aux_idx, 'UniformOutput', false);
end
realizations = cell2mat(reshape(realizations(:),1,1,cant_anns));

