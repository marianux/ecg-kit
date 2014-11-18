%% Remove baseline wandering with the median estimation method
%  This scripts exemplifies the baseline wander removal by
%  estimation/substraction and linear filtering.
%  For details refers to:
% 
%  Philip de Chazal; O'Dwyer, M.; Reilly, R.B - Automatic classification of
%  heartbeats using ECG morphology and heartbeat interval features. IEEE
%  TRANSACTIONS ON BIOMEDICAL ENGINEERING, VOL. 51, NO. 7, JULY 2004. Section
%  IV-A. ECG Filtering
%  
%  Sörnmo L, Laguna P. Bioelectrical Signal Processing in Cardiac
%  and Neurological Applications. Elsevier, 2005. ISBN
%  0-12-437552-9. Page 457.
%   
% Example
% 
%   [ECGblr, delayBLR]= BaselineWanderRemovalMedian( noisyECG, SamplingFreq)
% 
%   where:
%     *noisyECG is a the signal matrix, signals in columns. 
%     *SamplingFreq is the sampling frequency.
% 
%   and:
%     *ECGblr is the clean signal 
%     *delayBLR is the delay introduced by the method.
% 
% See also BaselineWanderRemovalSplines
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Birthdate: 01/01/2006
% Last update: 18/10/2014
% Copyright 2008-2014
% 
function [ECGblr, delayBLR]= BaselineWanderRemovalMedian( noisyECG, SamplingFreq)

WinSize = round(0.2*SamplingFreq);
delayBLR = round((WinSize-1)/2);

%Baseline estimation.
BaselineEstimation = MedianFilt(noisyECG, WinSize );

WinSize = round(0.6*SamplingFreq);
delayBLR = delayBLR + round((WinSize-1)/2);

BaselineEstimation = MedianFilt(BaselineEstimation, WinSize );

%Baseline removal
% ECGblr = zeros(size(noisyECG));
% 
% AuxRange = delayBLR:delayBLR+size(BaselineEstimation,1)-1;
% ECGblr(AuxRange,:) = noisyECG(AuxRange,:) - BaselineEstimation;

ECGblr = noisyECG - BaselineEstimation;



