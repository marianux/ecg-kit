function Output = PrctileFilt(Input, WinSize, Percentiles)

if(nargin < 2 || isempty(WinSize) )
    WinSize = 31;
end

MidPoint = ceil(WinSize/2);

aux_seq = 1:size(Input,1);
laux_seq = length(aux_seq);

each_sample_idx = arrayfun(@(a)(max(1,a):min(laux_seq,a + WinSize-1)), aux_seq - MidPoint+1, ...
                                                'UniformOutput', false);
Output = cellfun(@(a)(rowvec(prctile(Input(a,:), Percentiles, 1))), colvec(each_sample_idx), 'UniformOutput', false);

Output = cell2mat(Output);
