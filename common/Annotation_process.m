%% Convert heartbeat type of annotation from valid ECG formats to EC57 AAMI
% An alternative interface to AnnotationFilterConvert.
% 
% Example
% 
%   [QRS_locations, true_labels ] = Annotation_process(ann, recording_format, labeling_format)
% 
%   where:
%     *ann is a annotation structure (ECG_annotations prop. in ECGwrapper
%        object). Mandatory fields for this function are: time and anntyp. 
%     *recording_format is a string with a valid recording format read by
%        read_ecg (cKnownFormats).
%     *labeling_format is a string to map the annotations, possible values
%        are AAMI (N,S,V,F,Q) or AAMI2 (N,S,V&F,Q).
% 
% See also AnnotationFilterConvert
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Birthdate: 01/01/2012
% Last update: 18/10/2014
% Copyright 2008-2014
function [QRS_locations, true_labels ] = Annotation_process(ann, recording_format, labeling_format)

if( nargin < 3 || isempty(labeling_format) )
    %default AAMI2
    labeling_format = 'AAMI2';
end

%In this case we need the true labels labeled from an expert.
if( isfield( ann, 'anntyp'))
    % Annotations filtering and conversion
    ann_filtered = AnnotationFilterConvert(ann, recording_format, labeling_format);

    %debug and performance asssesment only 
    true_labels = ann_filtered.anntyp;
    QRS_locations = ann_filtered.time;
    
else
    
    QRS_locations = ann.time;
    true_labels = [];
    
end
