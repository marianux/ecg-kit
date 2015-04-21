%% (Internal) Mean/Median filtering
%   
%   Output = MedianFilt(Input, WinSize, bRobust)
% 
% Arguments:
% 
%      + Input: the signal
% 
%      + WinSize: The size of the window in samples
% 
%      + bRobust: use mean (false) or median (true)
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
%     BaselineEstimation = MedianFilt(noisyECG, WinSize );
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
function Output = MedianFilt(Input, WinSize, bRobust)

if(nargin < 2 || isempty(WinSize) )
    WinSize = 31;
end

if(nargin < 3 || isempty(bRobust) )
    bRobust = true;
end

if(bRobust)
    mean_ptr_func = @nanmean;
else
    mean_ptr_func = @nanmean;
end

MidPoint = ceil(WinSize/2);

aux_seq = 1:size(Input,1);
laux_seq = length(aux_seq);

each_sample_idx = arrayfun(@(a)(max(1,a):min(laux_seq,a + WinSize-1)), aux_seq - MidPoint+1, ...
                                                'UniformOutput', false);
Output = cellfun(@(a)(mean_ptr_func(Input(a,:),1)), colvec(each_sample_idx), 'UniformOutput', false);

Output = cell2mat(Output);


%old version
% startSample = MidPoint;
% endSample = size(Input,1) - fix(WinSize/2);
% Output = zeros(endSample-startSample+1, size(Input,2));
% 
% iCount = 1;
% 
% for iSampleIndex = startSample:endSample
%     
%     startRange = iSampleIndex-MidPoint+1;
%     iRange = startRange:startRange+WinSize-1;
%     
%     Output(iCount,:) = median( Input(iRange,:) );
%     
%     iCount = iCount + 1;
%     
% end
