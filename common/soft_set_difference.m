function [sft_set_difference idx]= soft_set_difference( idx1, idx2, win_size )

% Find the set difference between index sequence idx1-2 within a windows
% win_size.

sft_intersect = soft_intersect( idx1, idx2, win_size );

[sft_set_difference idx] = setdiff(idx1, sft_intersect);
