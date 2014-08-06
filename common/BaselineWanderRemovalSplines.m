function [ ECG BaselineWander ] = BaselineWanderRemovalSplines( ECG, QRS_locations, sampling_rate)

% Description:
% Performs baseline wander removal with cubic splines method. Estimates the
% baseline wander in the PQ silence segment and then substract it from ECG.
% 
% Arguments:
%   + ECG: signal to be cleaned.
%   + sampling_rate: sampling rate of the ECG.
% 
% Output:
%   + ECG: clean ECG.
% 
% References:
% 
% S�rnmo L, Laguna P. Bioelectrical Signal Processing in Cardiac
% and Neurological Applications. Elsevier, 2005. ISBN
% 0-12-437552-9. Page 457.
% 
% Limits and Known bugs:
% 
% Author: Mariano Llamedo Soria (llamedom at {electron.frba.utn.edu.ar; unizar.es}
% Birthdate  : 23/8/2011
% Last update: 23/8/2011

[ECG_size ECG_leads] = size(ECG);

QRS_locations = QRS_locations(QRS_locations > round(0.2*sampling_rate) & QRS_locations <= ECG_size );

cant_QRS = length(QRS_locations);

aux_idx = arrayfun(@(a)( QRS_locations(a) - round(0.2*sampling_rate): ...
                         min(ECG_size, QRS_locations(a) - round(0.1*sampling_rate))) , ...
                   1:cant_QRS, 'UniformOutput', false);

PQ_estimations = cell2mat(cellfun(@(a)(mean(ECG(a,:),1)), colvec(aux_idx), 'UniformOutput', false));
PQ_sample = cell2mat(cellfun(@(a)(mean(a)), colvec(aux_idx), 'UniformOutput', false));

% as spline is very memory-inneficient, we resample the PQ guide to ~25 Hz
sampling_rate_low = 25; % Hz
[ interp_k, deci_k ] = rat( sampling_rate_low / sampling_rate );
sampling_ratio = interp_k / deci_k;
ECG_size_decim = floor(ECG_size*sampling_ratio);

PQ_sample = round(PQ_sample * sampling_ratio );

eq_samples = find(diff(PQ_sample) == 0);

PQ_sample(eq_samples) = [];
PQ_estimations(eq_samples,:) = [];

if(PQ_sample(end) == ECG_size_decim)
    PQ_sample = PQ_sample(1:end-1);
    PQ_estimations = PQ_estimations(1:end-1,:);
end
if(PQ_sample(1) == 1)
    PQ_sample = PQ_sample(2:end);
    PQ_estimations = PQ_estimations(2:end,:);
end

% Spline diverges in many cases and could not fix it.
% BaselineWander = spline( [ 1       ; PQ_sample      ; ECG_size_decim  ], ...
%                 [ ECG(1,:); PQ_estimations ; ECG(end,:)]', ...
%                1:ECG_size_decim )';
BaselineWander = interp1(   [ 1       ; PQ_sample      ; ECG_size_decim  ], ... 
                            [ ECG(1,:); PQ_estimations ; ECG(end,:)], ...
                            1:ECG_size_decim, 'pchip' );

% substract the upsampled signal
ECG_size_new = ceil(ECG_size_decim/sampling_ratio);
last_samples = ECG(ECG_size_new,:);
ECG(1:ECG_size_new,:) = ECG(1:ECG_size_new,:) - resample(BaselineWander, deci_k, interp_k);
ECG(ECG_size_new+1:end,:) = bsxfun( @minus, ECG(ECG_size_new+1:end,:), (last_samples - ECG(ECG_size_new,:)));
