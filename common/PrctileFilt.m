%% (Internal) Arbitrary percentile filtering
%   
%   Output = PrctileFilt(Input, WinSize, Percentiles)
% 
% Arguments:
% 
%      + Input: the signal
% 
%      + WinSize: The size of the window in samples
% 
%      + Percentiles: the percentile to calculate, from 0 to 100, being
%      them the output signal of the filter.
% 
% Output:
% 
%      + Output: filtered output
% 
% Example:
% 
%     WinSize = round(0.2*SamplingFreq);
% 
%     %Baseline estimation.
%     BaselineEstimation = PrctileFilt(noisyECG, WinSize, );
% 
% 
% See also BaselineWanderRemovalMedian
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
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
