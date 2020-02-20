%% (Internal) Create new QRS detections based on other lead/algorithms detections
%   
%   artificial_annotations = combine_anns(time_serie, estimated_labs, header)
% 
% Arguments:
% 
%      + time_serie: all the QRS detections 
% 
%      + estimated_labs: estimated labels (TP, FP, FN) of each detection
% 
%      + header: ECG header struct.
% 
% Output:
% 
%      + artificial_annotations: all the QRS detections plus the new ones
%      created.
% 
% Example:
% 
%         aux_idx = best_detections_idx(1:min(10, length(best_detections_idx)));
%         artificial_annotations = combine_anns(all_annotations(aux_idx), estimated_labs(aux_idx), ECG_struct.header);
% 
% See also QRScorrector
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
function artificial_annotations = combine_anns(time_serie, estimated_labs, header)

    % cantidad de leads artificiales que generar�
    min_artificial_leads = 3;
    win_size = 20; % seconds
    win_size = round(win_size*header.freq); % milliseconds

%     lreferences = length(time_serie);
    for ii = 1:min_artificial_leads
        artificial_annotations(ii).time = [];
    end

    start_sample = min(cell2mat(cellfun(@(a)(min(a)),time_serie, 'UniformOutput', false)));
    end_sample = max(cell2mat(cellfun(@(a)(max(a)),time_serie, 'UniformOutput', false))) + 1;
    
%     aux_seq = (start_sample+win_size):round(win_size/2):end_sample;
    aux_seq = min(start_sample+win_size, end_sample):win_size:end_sample;
    
    if( isempty(aux_seq) )
        str_aux = disp_string_framed(0, 'Recording too short for combining detections'  );
        fprintf(2, '\n%s\n', str_aux);
        return
    end
    
    if(aux_seq(end) ~= end_sample )
        aux_seq = [aux_seq (aux_seq(end)+win_size) ];
    end
    
    laux_seq = length(aux_seq);
    
    aux_idx = arrayfun(@(a)( ...
                                cellfun( @(b)( ...
                                                findStartEnd(  b >= (a - win_size) & b < a ) ...
                                            ), time_serie, 'UniformOutput', false  )... 
                            ), aux_seq, 'UniformOutput', false);
    
    aux_q = cellfun( @(a)( cell2mat(cellfun( @(b, c)( ...
                                    calc_q_val( b, c ) ...
                                ), estimated_labs, a, 'UniformOutput', false  ))... 
                ), aux_idx, 'UniformOutput', false  );

    [~, aux_q_idx] = cellfun( @(a)( sort(a, 'descend') ), aux_q, 'UniformOutput', false );

    % TODO: hacer una estrategia que divida los segmentos de  aux_seq hasta
    % que no haya problemas de transiciones de los segmentos solapados.
%     while()
%         aux_idx = arrayfun( @(ii)( cell2mat(arrayfun( @(jj)( find_disagreements(aux_idx, aux_q_idx, ii, jj ) ), 1:laux_seq , 'UniformOutput', false)) ), 1:3, 'UniformOutput', false );
%     end
        
    aux_val = arrayfun( @(ii)( cell2mat(cellfun( @(a, q_idx)( build_combined_series(time_serie, a, q_idx, ii ) ), aux_idx, aux_q_idx, 'UniformOutput', false)) ), 1:min_artificial_leads, 'UniformOutput', false );
    
    cant_artificial_leads = min(min_artificial_leads, sum(cellfun(@(a)(~isempty(a)), aux_val )) );
    
    if( cant_artificial_leads > 0 ) 
        
        for ii = 1:cant_artificial_leads
            % avoid annotations very close each other.
            aux_time_serie = aux_val{ii};
            aux_time_serie(find( diff(sort(aux_time_serie)) <= round(0.15 * header.freq) ) +1) = [];
            if( ~isempty(aux_time_serie) )
                artificial_annotations(ii).time = colvec(aux_time_serie);
            end
        end
        
        if( cant_artificial_leads < min_artificial_leads )
            for ii = cant_artificial_leads+1:min_artificial_leads
                artificial_annotations(ii).time = artificial_annotations(cant_artificial_leads).time;
            end
        end
        
    end
    
function start_end_aux = findStartEnd( bAux )

    start_aux = find(  bAux, 1, 'first' );
    end_aux = find(  bAux, 1, 'last' );
    start_end_aux = [start_aux end_aux];


function this_q = calc_q_val(this_labs, strt_end)
%         win_size in milliseconds

    if( isempty(strt_end) || isempty(this_labs) )
        this_q = [];
    else
        this_labs = this_labs(strt_end(1):strt_end(2));

        this_se = sum(this_labs == 3) / sum(this_labs == 3 | this_labs == 1) ;
        this_pp = sum(this_labs == 3) / sum(this_labs == 3 | this_labs == 2) ;

        % paper metric
        %this_q = (2*this_se + this_pp)/3;
        % F1 score
        this_q = 2/(1/this_se + 1/this_pp);
    end
    
function new_str_end = find_disagreements(str_end, q_idx, ii, jj)
      
    this_q_idx = q_idx{jj};
    next_q_idx = q_idx{jj+1};
    
    if( isempty(this_q_idx) || ii > length(q_idx) )
        new_str_end = [];
    else
        if( q_idx(ii) == q_idx(ii+1) )
            new_str_end = [];
        else
            
        end
    end
    

function combined_series = build_combined_series(ts, str_end, q_idx, ii)
      
    if( isempty(q_idx) || ii > length(q_idx) )
        combined_series = [];
    else
        q_idx = q_idx(ii);
        aux_idx = str_end{q_idx};
        aux_val = ts{q_idx};
        if( isempty(aux_idx) || isempty(aux_val)  )
            combined_series = [];
        else
            combined_series = rowvec(aux_val(aux_idx(1):aux_idx(2)));
        end
    end
    
        
