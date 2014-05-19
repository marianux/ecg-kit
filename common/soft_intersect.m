function [sft_intersect1 idx1 idx2 ] = soft_intersect( val1, val2, win_size )

% Find the intersection between index sequence val1-2 within a windows
% win_size.
idx1 = [];
idx2 = [];

sft_intersect1 = unique(cell2mat(arrayfun(@(a)( colvec( val1(val1 >= (a - win_size) & val1 <= (a + win_size) ) ) ), colvec(val2), 'UniformOutput', false)));

if( nargout > 1)
    [ ~, idx1] = intersect(val1, sft_intersect1);
end

if( nargout > 2)
    sft_intersect2 = unique(cell2mat(arrayfun(@(a)( colvec( val2(val2 >= (a - win_size) & val2 <= (a + win_size) ) ) ), colvec(sft_intersect1), 'UniformOutput', false)));
    [ ~, idx2] = intersect(val2, sft_intersect2);
end
