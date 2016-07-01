%% (Internal) Pack a signal in a realization array, based on a sync. signal and a time window around.
%
%  [realizations, aux_idx, seq_idx] = pack_signal(signal, anns_refined, realization_limit, substract_mean)
% 
% Description:
% Packs several realizations synchronized at each sample in anns_refined
% (QRS detections), in a window from anns_refined(i) - realization_limit(1)
% to anns_refined(i) + realization_limit(2)
% 
% 
% Arguments:
% 
%      + signal: the signal
% 
%      + anns_refined: the synch sample to stack realizations
%             
%      + realization_limit: the limits of each realization, from
%      anns_refined(i) - realization_limit(1) 
%      to 
%      anns_refined(i) + realization_limit(2)  
% 
%      + substract_mean: Boolean to substract mean from each realization.
%             
% Output:
% 
%      + realizations : realizations is a signal pack where dimensions obtained
%      with size indicates:
%           [ ensemble_size nsig cant_anns ] = size(realizations);
% 
%      + aux_idx: the signal samples corresponding to realizations
% 
%      + seq_idx: the anns indexes corresponding to realizations
% 
% Example:
% 
% 
% See also ECGtask_ECG_delineation
% 
% Author: Mariano Llamedo Soria (llamedom at frba.utn.edu.ar)
% Version: 0.1 beta
% Birthdate  : 30/7/2014
% Last update: 30/7/2014
% Copyright 2008-2015
% 
function [realizations, aux_idx, seq_idx] = pack_signal(signal, anns_refined, realization_limit, substract_mean)

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

