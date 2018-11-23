%% (Internal) Find modulus maxima in a signal
%   
% Function which returns the indexes of vector x in which there are
% local modulus maxima or minima whose modulus is greater than 
% threshold.
% if signo is 0, it doesn't matter, if signo is +1 or -1, it only searchs
% for modulus maxima positive or negative
% t_restriction: time restriction. Doesnt allow modmax's closer than
% n_greater. In order to restrict the amount of detections, just keep the
% n_greater max_mod found. If more than one modmax, the bigger wins. 
% pb a progress bar object to track evolution of the algorithm in big
% signals
% 
%   [indexes max_mod] = modmax(x, first_samp, threshold, signo, t_restriction, n_greater)
% 
% Arguments:
% 
%      + x: the signal
% 
%      + first_samp: analyze signal from first_samp sample
%             
%      + threshold: an amplitude threshold to consider maxima
% 
%      + signo: sign of the maxima
%             
%      + t_restriction: time restriction between adjacent maximums
% 
%      + n_greater: return only the n_greater maxima
%             
% Output:
% 
%      + indexes : the indexes of the modulus maxima
% 
%      + max_mod: the values of the modulus maxima
% 
% Example:
% 
% Author: Juan Pablo Mart�nez, Rute Almeida
% adapted by Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
function [indexes max_mod] = modmax(x, first_samp, threshold, signo, t_restriction, n_greater, arg_pb)
lx = size(x,1);
indexes = [];
max_mod = [];

if( nargin < 2 || isempty(first_samp) )
    first_samp = [2 lx];
else
    if( length(first_samp) < 2 )
        first_samp = [max(2,first_samp) lx];
    else
        first_samp(1) = max(2,first_samp(1));
        first_samp(2) = min(lx,first_samp(2));
    end
end

if( nargin < 3 || isempty(threshold) )
    threshold = 0;
end

if( nargin < 4 || isempty(signo) )
    signo = 0;
end

bWithPB = false;
if( lx > 10e5 )
    bWithPB = true;
    if( nargin < 7 || isempty(arg_pb) )
        pb = progress_bar('Modmax function');
    else
        pb = arg_pb;
    end
end


if( lx > first_samp(1) )

    s = sign(x);
    x = abs(x);

    sample_curr_idx = first_samp(1):first_samp(2)-1;
    sample_prev_idx = (first_samp(1)-1):first_samp(2)-2;
    sample_next_idx = (first_samp(1)+1):first_samp(2);
    
    localmax =    ( x(sample_curr_idx,:)  >= x(sample_prev_idx,:) ) ...
                & ( x(sample_curr_idx,:)  >  x(sample_next_idx,:) ...
                &   x( sample_curr_idx,:) >= threshold) ...
                & ( s(sample_curr_idx,:)*signo >=0 );   % if 0,it doesnt matter

    iAux = false(size(x));
    iAux(sample_curr_idx,:) = localmax;
    indexes = find(iAux);
    max_mod = x(indexes) .* s(indexes);

else
    return
end

if( nargin < 5 || isempty(t_restriction) )
    t_restriction = 0;
end

if( t_restriction > 0 )
    
    aux_idx = [];
    ii = 1;
    bRefine = false; 

    if( bWithPB )
        pb.reset();
    end

    lindexes = length(indexes);
%     if( bWithPB )
%         pb.Loops2Do = lindexes;
%     end

    [~, aux_sorted_mod_idx] = sort(x(indexes), 'descend');

    while( ii < lindexes )

% too expensive computationally, check it !!
%         if( bWithPB )
%             pb.start_loop();
%         end
        
        if( ~isnan(indexes(aux_sorted_mod_idx(ii))) )

            indexes_inside_idx = find( indexes >= (indexes(aux_sorted_mod_idx(ii)) - t_restriction) & indexes <= (indexes(aux_sorted_mod_idx(ii))+t_restriction) );
            indexes_inside_idx(indexes_inside_idx == aux_sorted_mod_idx(ii)) = [];

            if( ~isempty(indexes_inside_idx) )
                indexes(indexes_inside_idx) = nan;
            end
        end

        ii = ii+1;

%         if( bWithPB )
%             pb.LoopsDone = ii;
%             pb.checkpoint('');
% 
%             pb.end_loop();
%         end
        
    end

    indexes(isnan(indexes)) = [];
        
    max_mod = x(indexes) .* s(indexes);
    
end

% if( bWithPB )
%     clear pb
% end


lindexes = length(indexes);

if( nargin < 6 || isempty(n_greater) )
    n_greater = lindexes;
end

if( n_greater < lindexes )
    
    [~, aux_idx] = sort(abs(max_mod), 'descend');
    
    indexes = indexes(aux_idx(1:n_greater));
    max_mod = max_mod(aux_idx(1:n_greater));

    [indexes, aux_idx] = sort(indexes);
    max_mod = max_mod(aux_idx);
    
end

