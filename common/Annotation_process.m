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
