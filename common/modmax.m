function [indexes max_mod] = modmax(x, first_samp, threshold, signo, t_restriction, n_greater, pb)
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

if( nargin < 7 || isempty(pb) )
    pb = progress_bar('Modmax function');
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

end

if( nargin < 5 || isempty(t_restriction) )
    t_restriction = 0;
end

if( t_restriction > 0 )
    
    bContinue = true; 
    while( bContinue )
        
        
        aux_idx = [];
        ii = 1;
        bRefine = false; 
        
        pb.reset();
        
        lindexes = length(indexes);
        pb.Loops2Do = lindexes;
        
        while( ii < lindexes )
            
            pb.start_loop();

            indexes_inside_idx = find( indexes >= indexes(ii) & indexes < (indexes(ii)+t_restriction) );
            if( length(indexes_inside_idx) == 1 )
                aux_idx = [ aux_idx; indexes(ii)];
                ii = ii+1;
            else            
                bRefine = true; 
                max_idx = max_index( x(indexes(indexes_inside_idx)) );
                aux_idx = [ aux_idx; indexes(indexes_inside_idx(max_idx))];
                ii = indexes_inside_idx(end)+1;
            end

            pb.LoopsDone = ii;
            pb.checkpoint('');
            
            pb.end_loop();
            
        end

        indexes = aux_idx;
        
        if( ~bRefine )
            bContinue = false; 
        end
        
    end
    
    max_mod = x(indexes) .* s(indexes);
    
end


if( nargin < 6 || isempty(n_greater) )
    n_greater = realmax;
end

lindexes = length(indexes);

if( n_greater < lindexes )
    
    [~, aux_idx] = sort(abs(max_mod), 'descend');
    
    indexes = indexes(aux_idx(1:n_greater));
    max_mod = max_mod(aux_idx(1:n_greater));

    [indexes, aux_idx] = sort(indexes);
    max_mod = max_mod(aux_idx);
    
end

