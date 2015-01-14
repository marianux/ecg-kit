%% Fetch the best QRS detections from an ECGtask_QRSdetections object
% Get the best QRS detections according to the ratios quality metric
% provided by the QRS detection task.
%   
% Example
% 
%   [QRS_detections, str_delineation_chosen] = GetBestQRSdetections(payload)
% 
% Arguments:
%      +payload: [struct] REQUIRED
%           
%           An struct from loading the resulting mat file of a QRS
%           detection task.
% 
% Output:
%     + QRS_detections: the detections indexes
% 
%     + str_delineation_chosen: The name or id of the QRS detection indexes.
% 
% See also ECGwrapper
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 28/11/2014
% Birthdate  : 28/11/2014
% Copyright 2008-2014
function [QRS_detections, str_delineation_chosen] = GetBestQRSdetections(payload)

QRS_detections = [];
str_delineation_chosen = [];

aux_fn = fieldnames(payload);
% prefer manually reviewed annotations corrected 
aux_idx = find(cell2mat( cellfun(@(a)(~isempty(strfind(a, 'corrected_'))), aux_fn, 'UniformOutput', false)));

if( isempty(aux_idx) )

    if( isfield(payload, 'series_quality') )
        % choose the best ranked automatic detection
        [~, aux_idx] = sort(payload.series_quality.ratios, 'descend' );
        QRS_detections = payload.(payload.series_quality.AnnNames{aux_idx(1),1} ).(payload.series_quality.AnnNames{aux_idx(1),2});
        str_delineation_chosen = payload.series_quality.AnnNames{aux_idx(1),1};
    else
        % choose the first of wavedet
        aux_idx = find(cell2mat( cellfun(@(a)(~isempty(strfind(a, 'wavedet_'))), aux_fn, 'UniformOutput', false)));
        if( isempty(aux_idx) )
            disp_string_framed(2, 'Could not identify QRS detections in the input payload.');
            return
        else
            QRS_detections = payload.(aux_fn{aux_idx(1)}).time;
            str_delineation_chosen = aux_fn{aux_idx(1)};
        end
    end
else
    % choose the first manually audited
    QRS_detections = payload.(aux_fn{aux_idx(1)}).time;
    str_delineation_chosen = aux_fn{aux_idx(1)};
end
