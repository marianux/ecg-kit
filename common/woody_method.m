%% (Internal) Woody algorithm for heartbeat allignment
%
%     [ best_esemble_avg, 
%       min_anns_refined, 
%       all_anns_matched_with_BEA] = woody_method(signal, 
%                           annotations, realization_limit, 
%                           outliers_proportion, bRobust )
% 
% Arguments:
% 
%   + signal:
% 
%   + annotations: heartbeat sample location.
% 
%   + realization_limit: limits of the heartbeat with respect to
%   annotations(i) - realization_limit(1) to annotations(i) +
%   realization_limit(2).  
% 
%   + outliers_proportion: Proportions of true heartbeats respect to
%   outliers to assume, and discard from the template calculation. Default 0.9.  
% 
%   + bRobust: Calculate the template with median instead of mean.
%
%   + bSubstract_mean: Remove mean value from each realization.
% 
% Output:
% 
%   + best_esemble_avg: The best template found.
% 
%   + min_anns_refined: The refined location of the selected heartbeats,
%   without outliers. 
% 
%   + all_anns_matched_with_BEA: The refined location of the selected heartbeats,
%   with outliers. 
% 
% See also 
% 
% Author: Mariano Llamedo Soria (llamedom at frba.utn.edu.ar)
% Version: 0.1 beta
% Birthdate: 21/7/2010
% Last update: 20/02/2013
% Copyright 2008-2015
% 
function [best_esemble_avg, min_anns_refined, all_anns_matched_with_BEA] = woody_method(signal, annotations, realization_limit, outliers_proportion, bRobust, bSubstract_mean )

if( nargin < 6 || isempty(bSubstract_mean ) )
    bSubstract_mean = true;
end

if( nargin < 5 || isempty(bRobust ) )
    bRobust = false;
end

MAX_LOOPS = 20;

[nsamp, nsig] = size(signal);

cant_anns = length(annotations);
anns_refined = colvec(annotations);
discarded_anns_idx = [];
anns_aux_idx = colvec(1:cant_anns);
loop_count = 1;
crit = [];
crit(1) = inf;

if( nargin < 4 || isempty(outliers_proportion) )
    outliers_proportion = 0.9;
end

min_cant_anns = round(outliers_proportion * cant_anns);

ensemble_size = sum(realization_limit);

min_anns_refined = [];
min_crit = inf;


while( loop_count < MAX_LOOPS && crit(max(loop_count-1,1)) > 0.001 )

    [avg_pack, aux_idx, anns_idx ] = pack_signal(signal, anns_refined, realization_limit, bSubstract_mean);    

    if( bRobust )
        esemble_avg = flipud(squeeze(median(avg_pack,3)));
    else
        esemble_avg = flipud(squeeze(mean(avg_pack,3)));
    end

    if(loop_count == 1) 
        esemble_avg_initial = esemble_avg;
    end

    realizations = cellfun(@(a)( cell2mat(arrayfun( @(aa)( conv(signal(a,aa), esemble_avg(:,aa) ) ), 1:nsig, 'UniformOutput', false )) ), aux_idx, 'UniformOutput', false);
    
    cross_corr_idx = cellfun(@(a)( arrayfun( @(aa)( max_index( a(:,aa) ) ), 1:nsig)), realizations, 'UniformOutput', false);
    cross_corr_idx = cell2mat( cross_corr_idx(:) ) ;
    
    %just calculate the median change in time
%     error_sig = round(median( cross_corr_idx, 2 )) - ensemble_size ;
    
    %do a weighted avg favoring small changes in time
    cross_corr_idx = cross_corr_idx - ensemble_size;
    err_weight = pdf('norm',cross_corr_idx,0,5); % 5 samples of refinement
    error_sig = round( sum(cross_corr_idx .* err_weight,2) ./ (sum(err_weight,2)+1e-6) );

% % debugging lines
% plot_ecg_mosaic(shiftdim(avg_pack,1), [], [], [], [], [], 3);
% plot_ecg_mosaic(flipud(esemble_avg), [], [], [], [], [], 4);
% realizations_discarded = pack_signal(signal, annotations(discarded_anns_idx), realization_limit);    
% plot_ecg_mosaic(shiftdim(realizations_discarded,1), [], [], [], [], [], 5);
% [~, greater_change_idx] = sort(abs(error_sig), 'descend');
% realizations_greater_change = pack_signal(signal, anns_refined(greater_change_idx(1:10)), realization_limit);    
% plot_ecg_mosaic(shiftdim(realizations_greater_change,1), [], [], [], [], [], 6);
% figure(7)
% hist(realizations_out)
    
    % time correction
    anns_refined(anns_idx) = anns_refined(anns_idx) + error_sig;

    crit(loop_count) = sqrt(mean(error_sig.^2));

    if(crit(loop_count) < min_crit)

        min_crit = crit(loop_count);
        min_anns_refined = anns_refined;
        min_discarded_anns_idx = discarded_anns_idx;
        best_esemble_avg = flipud(esemble_avg);

    end

    if( cant_anns > min_cant_anns )
        %try to refine the realizations with the most similar examples
        aux_val = cell2mat(reshape(realizations, 1, 1, length(realizations))); 
        realizations_lims = prctile( aux_val, [2.5 97.5], 3);
        realizations_out = bsxfun( @lt, aux_val, realizations_lims(:,:,1)) | bsxfun( @ge, aux_val, realizations_lims(:,:,2));
        realizations_out = squeeze(sum(sum(realizations_out, 1),2));
        [~, realizations_idx] = sort( realizations_out, 'Descend' );

        %discard most outlying examples
        aux_idx = colvec(find(anns_refined < realization_limit(1) | anns_refined > (nsamp - realization_limit(2))));
        aux_idx = unique([aux_idx; realizations_idx(1:round(0.05*cant_anns) )]);
        discarded_anns_idx = [discarded_anns_idx; anns_aux_idx( aux_idx ) ];
        anns_refined(aux_idx ) = [];
        anns_aux_idx(aux_idx ) = [];
        cant_anns = length(anns_refined);
        
    else
        % watch refinements not move to the extreme of the segment
        aux_idx = colvec(find(anns_refined < realization_limit(1) | anns_refined > (nsamp - realization_limit(2))));
        discarded_anns_idx = [discarded_anns_idx; anns_aux_idx( aux_idx ) ];
        anns_refined(aux_idx ) = [];
        anns_aux_idx(aux_idx ) = [];
        cant_anns = length(anns_refined);
    end
    
    loop_count = loop_count + 1;
end

difference_between_pattern = cell2mat(arrayfun( @(aa)( conv(esemble_avg_initial(:,aa), best_esemble_avg(:,aa) )/(size(esemble_avg_initial,1)*std(best_esemble_avg(:,aa))* std(esemble_avg_initial(:,aa)) ) ), 1:nsig, 'UniformOutput', false ));

cross_corr_val = mean(difference_between_pattern(ensemble_size,:));

if( loop_count == MAX_LOOPS && cross_corr_val < 0.8 )
    
    if( usejava('desktop') )
        
        %warn humans.
        warning('woody_method_BAD_CONVERGENCE', 'Woody method reach the last iteration with a poor result. Check results manually.')

%     else
%         best_esemble_avg = [];
%         min_anns_refined = annotations;
%         all_anns_matched_with_BEA = annotations;
%         return
    end
    
end

if( nargout > 2 )

    all_anns_discarded = sort([ annotations(min_discarded_anns_idx)] );

    cant_anns = length(all_anns_discarded);
    seq_idx = 1:cant_anns;

    aux_idx = arrayfun(@(a)(all_anns_discarded(a) - realization_limit(1):all_anns_discarded(a) + realization_limit(2)), ...
                            seq_idx, 'UniformOutput', false);

    realizations = cellfun(@(a)( cell2mat(arrayfun( @(aa)( conv(signal(max(1,min(nsamp,a)),aa), flipud(best_esemble_avg(:,aa)) ) ), 1:nsig, 'UniformOutput', false )) ), aux_idx, 'UniformOutput', false);

    cross_corr_idx = cellfun(@(a)( arrayfun( @(aa)( max_index( a(:,aa) ) ), 1:nsig)), realizations, 'UniformOutput', false);
    cross_corr_idx = cell2mat( cross_corr_idx(:) ) ;
    error_sig = round(median( cross_corr_idx, 2 )) - ensemble_size ;

    all_anns_discarded = colvec(all_anns_discarded) + colvec(error_sig);

    all_anns_matched_with_BEA = sort([ min_anns_refined; all_anns_discarded] );

end


