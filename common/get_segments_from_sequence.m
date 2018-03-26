%% (Internal) Get the time samples of a (xx,yy) sequence below two restrictions (x_thr, y_thr)
%
% This function calculates the time samples (xx) of a two valued sequence
% (xx, yy) below two restrictions (x_thr, y_thr)
% 
%       selected_segments = get_segments_from_sequence(yy, xx, y_thr, x_thr )
% 
% Arguments:
% 
%       xx, yy: xx and yy values. 
% 
%       x_thr, y_thr: the thresholds to search for within the sequence.
% 
% Output:
% 
%       selected_segments: a (n x 2) matrix with the segments start-end
%                          where n is the amount of segments found in the
%                          sequence 
% 
% Example
% 
% See also RR_calculation, MedianFiltSequence
% 
% Author: Mariano Llamedo Soria (llamedom at frba.utn.edu.ar)
% Version: 0.1 beta
% Birthdate  : 9/5/2017
% Last update: 9/5/2017
% Copyright 2008-2017
% 
function selected_segments = get_segments_from_sequence(xx, yy, x_thr, y_thr )

    selected_segments = [];

    % sequence values below the thr
    bAux = yy <= y_thr;
    % indexes above the thr
    aux_idx = find(~bAux);
    
    aux_val = diff(xx(aux_idx));
    % starts of segments at least x_thr long that are below y_thr
    aux_idx3 = find(aux_val >= x_thr);
    below_idx = xx(bAux);
    above_idx = xx(~bAux);

% figure(3); plot(xx, yy ); xlims = xlim(); ylims = ylim(); ylims = ylims + [0.1 -0.1] * diff(ylims); hold on; plot(xlims, [y_thr y_thr] , '--r' ); hold off;

    for jj = 1:length(aux_idx3) 
        start_idx = xx(aux_idx(aux_idx3(jj))+1);
        end_idx = find( above_idx > start_idx, 1,'first');
        if isempty(end_idx)
            continue
        end
        end_idx = below_idx( find( below_idx < above_idx(end_idx) ,1, 'last') );
        
        if( ~isempty(start_idx) && ~isempty(end_idx) )
            selected_segments = [selected_segments; [start_idx end_idx]];
        end
    end
    
end

