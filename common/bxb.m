function [C, TP_ann_idx, TP_det_idx, FN_idx, FP_idx] = bxb( ann, test_qrs_idx, heasig, strDBformat)

if( nargin < 4 || isempty(strDBformat) )
    strDBformat = 'MIT';
end
    
AAMI_anntypes =  ('NSVFQ')';

%filtramos las anotaciones que son latidos

ann = AnnotationFilterConvert(ann, strDBformat, 'AAMI');
ann.anntyp = AAMI_anntypes(ann.anntyp);

[test_qrs_idx, test_qrs_idx2] = sort(test_qrs_idx);

%find matchs now
lRefBeats = length(ann.anntyp);
lDetectedBeats = length(test_qrs_idx);
Ref_idx = 1;
win_in_samples = round(0.15*heasig.freq);
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
