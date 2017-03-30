%% (Internal) Compares two heartbeat series and produce the confusion matrix as result
%   
%   [C, TP_ann_idx, TP_det_idx, FN_idx, FP_idx] = bxb( ann, test_qrs_idx, heasig, strDBformat)
% 
% Arguments:
% 
%      + ann: the string
% 
%      + test_qrs_idx: size of the target string.
% 
%      + heasig: "left" "right" "center"
% 
%      + strDBformat: the string
% 
% Output:
% 
%      + C: the confusion matrix.
% 
%      + TP_ann_idx, TP_det_idx, FN_idx, FP_idx: the indexes of the true
%      positive (TP), false negative (FN) and FP within test_qrs_idx
% 
% Example:
% 
% See also ECGtask_QRS_detection
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
function [C, TP_ann_idx, TP_det_idx, FN_idx, FP_idx] = bxb( ann, test_qrs_idx, heasig, strDBformat, bFilterAnnotations)

% constants
win_in_samples = round(0.15*heasig.freq);


if( nargin < 5 || isempty(bFilterAnnotations) )
    bFilterAnnotations = false;
end

if( nargin < 4 || isempty(strDBformat) )
    strDBformat = 'MIT';
end

if( isstruct(test_qrs_idx) && isfield(ann, 'time') )
    test_qrs_idx = test_qrs_idx.time;
end

%filtramos las anotaciones que son latidos
if( bFilterAnnotations && isfield(ann, 'anntyp') )
    ann = AnnotationFilterConvert(ann, strDBformat, 'AAMI');
end


% AHA database special case: just annotated a small part of each recording.
if( strcmpi(strDBformat, 'AHA') )
    test_qrs_idx = test_qrs_idx(test_qrs_idx >= (ann.time(1) - win_in_samples) & test_qrs_idx <= (ann.time(end) + win_in_samples));
end


[test_qrs_idx, test_qrs_idx2] = sort(test_qrs_idx);

%find matchs now
lRefBeats = length(ann.time);
lDetectedBeats = length(test_qrs_idx);
Ref_idx = 1;
TP_idx = [];
FN_idx = [];

while( Ref_idx <= lRefBeats )
    
    beat_distance = abs( ann.time(Ref_idx) - test_qrs_idx);
    close_beats_idx = find( beat_distance <= win_in_samples );
    
    if( ~isempty(close_beats_idx) )
        
        lclose_beats = length(close_beats_idx);
        
        if( lclose_beats == 1  )
            %matched beat
            TP_idx = [TP_idx; Ref_idx test_qrs_idx2(close_beats_idx)];
        
        elseif( lclose_beats > 1  )
            %several beats close to the true beat.
            [~, closest_beat_idx] = sort( beat_distance(close_beats_idx) );
            
            if( beat_distance(close_beats_idx(closest_beat_idx(1))) ~= beat_distance(close_beats_idx(closest_beat_idx(2))) )
                %choose the closest
                TP_idx = [TP_idx; Ref_idx test_qrs_idx2(close_beats_idx(closest_beat_idx(1))) ];
            else
                %choose the first happened match 
                if( test_qrs_idx(close_beats_idx(closest_beat_idx(1))) <= test_qrs_idx(close_beats_idx(closest_beat_idx(1))) )
                    TP_idx = [TP_idx; Ref_idx test_qrs_idx2(close_beats_idx(closest_beat_idx(1))) ];
                else
                    TP_idx = [TP_idx; Ref_idx test_qrs_idx2(close_beats_idx(closest_beat_idx(2))) ];
                end
            end
            
        end
    else
        %beat missed
        FN_idx = [FN_idx; Ref_idx];
    end
    
    Ref_idx = Ref_idx + 1;
    
end

%build the beat index references.
if( isempty(TP_idx) )
    TP_ann_idx = [];
    TP_det_idx = [];
else
    TP_ann_idx = TP_idx(:,1);
    TP_det_idx = TP_idx(:,2);
end

FP_idx = setdiff(1:lDetectedBeats, TP_det_idx);

%build the confusion matrix.
C = [ length(TP_ann_idx), length(FN_idx); length(FP_idx), 0 ];

if(nargout == 0)
        %pretty print performance
        
end
