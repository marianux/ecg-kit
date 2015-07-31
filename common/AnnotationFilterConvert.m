%% Convert heartbeat type of annotation from valid ECG formats to EC57 AAMI
% Map heartbeat annotations in the anntyp field of the ann structure to
% valid EC57 AAMI classes. Annotations that are not heartbeats are
% filtered out.
% 
% Example
% 
%   ann = AnnotationFilterConvert(ann, recording_format, labeling_format)
% 
%   where:
%     *ann is a annotation structure (ECG_annotations prop. in ECGwrapper
%        object). Mandatory fields for this function are: time and anntyp. 
%     *recording_format is a string with a valid recording format read by
%        read_ecg (cKnownFormats).
%     *labeling_format is a string to map the annotations, possible values
%        are AAMI (N,S,V,F,Q) or AAMI2 (N,S,V&F,Q).
% 
% See also ECGwrapper, read_ECG
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Birthdate: 01/01/2012
% Last update: 18/10/2014
% Copyright 2008-2015% Version: 0.1 beta
% Birthdate: 01/01/2012
% Last update: 18/10/2014
% Copyright 2008-2015
function ann = AnnotationFilterConvert(ann, recording_format, labeling_format)

cLabelingFormats = {'AAMI', 'AAMI2'};
LabelingIndexes = [2; 3];

if( nargin < 1 || isempty(ann) )
    ann = [];
    return
end

if( isfield(ann,'time') && isempty(ann.time) )
    return
end

if( nargin > 2 && ~isempty(labeling_format) )
    LabConversion_idx = LabelingIndexes(find(strcmpi(labeling_format, cLabelingFormats)));
else
    %default AAMI2
    LabConversion_idx = LabelingIndexes(2);
end



if( strcmpi(recording_format, 'AHA') )

%%%%%%%%%%%%%%%%%
%%% AHA FORMAT %%%
%%%%%%%%%%%%%%%%%%%


    % Esta matriz es OBSOLETA. Fue reemplazada por
    % iLabListsTranslations. La deje para documentacion de las
    % clases nada mas.
    cTranslationMatrix = { ...
        'Ventricular Escape'                                'V Escape'                        9   'E'  'E' ; ...
        'Fusion Beat'                                       'Fusion Beat'                     10  'F'  'F' ; ...
        'Beat of Non-Ventricular Origin'                    'Beat of Non-Ventricular Origin'  1   'N'  'N' ; ...
        'Paced Beat'                                        'Paced Beat'                      4   'P'  'P' ; ...
        'Questionable Beat - Indeterminate Origin'          'Questionable Beat'               5   'Q'  'Q' ; ...
        'R-on-T Beat'                                       'R-on-T Beat'                     8   'R'  'R' ; ...
        'Unreadable'                                        'Unreadable'                      5   'U'  'U' ; ...
        'Premature Ventricular Contraction'                 'Ventricular'                     8   'V'  'V' ; ...
        'Beginning of Ventricular Fibrillation or Flutter'  'Start Flutter o Fib'             9   '['  '[' ; ...
        'End of Ventricular Fibrillation or Flutter'        'End Flutter o Fib'               9   ']'  ']' ; ...
         'fv c/ NaN'                                        'fv c/ NaN'                       9   '�'  '�' ; ...
                        };

    cValidAnnotations =          cTranslationMatrix(:,5);
    cAAMItranslation =           cTranslationMatrix(:,4);
    cAnnTranslation   =          cTranslationMatrix(:,2);
    cValidAnnotationsString =    cTranslationMatrix(:,1);

    cLablists = { 'AHA' 'AAMI' 'AAMI2' 'Rrhythm' 'Morph '};

    cLabListsLabels = { ...
                            cValidAnnotations;
                            {'Normal' 'Supraventricular' 'Ventricular' 'Fusion' 'Unclass' 'NoAAMI'};
                            {'Normal' 'Supraventricular' 'Ventricular' 'Unclass' 'NoAAMI'};
                            {'Prematuro' 'Normal' 'Escape' 'Unclass'};
                            {'Normal' 'Aberrada' 'Unclass'} };


    iLabListsTranslations = [ ...
        ...% AHA labelling:    'E' 'F' 'N' 'P' 'Q' 'R' 'U' 'V' '[' ']' 'ñ'
                                1:size(cTranslationMatrix,1);
                                3   4   1   6   5   3   6   3   6   6   6;
                                3   3   1   5   4   3   5   3   5   5   5;
                                3   2   2   4   4   1   4   1   4   4   4;
                                2   1   1   3   3   2   3   2   3   3   3;
                            ];
                        
elseif( strcmpi(recording_format, 'HES') )

%%%%%%%%%%%%%%%%%%
%%% HES FORMAT %%%
%%%%%%%%%%%%%%%%%%


    % Esta matriz es OBSOLETA. Fue reemplazada por
    % iLabListsTranslations. La deje para documentacion de las
    % clases nada mas.
    cTranslationMatrix = { ...
                        'NO_TYPING'          'NO_TYPING'            5   'Q'  0   ; ...
                        'DOMINANT_TYPE_0'    'DOMINANT_TYPE_0'      1   'N'  1   ; ...
                        'DOMINANT_TYPE_1'    'DOMINANT_TYPE_1'      1   'N'  2   ; ...
                        'ABERRANT'           'ABERRANT'             2   'S'  4   ; ...
                        'ARTEFAKT'           'ARTEFAKT'             5   'Q'  5   ; ...
                        'VES'                'VES'                  3   'V'  60  ; ...
                        'SVES'               'SVES'                 2   'S'  70  ; ...
                        };

    cValidAnnotations =          cell2mat(cTranslationMatrix(:,5));
    cAAMItranslation =           cTranslationMatrix(:,4);
    cAnnTranslation   =          cTranslationMatrix(:,2);
    cValidAnnotationsString =    cTranslationMatrix(:,1);

    cLablists = { 'HES' 'AAMI' 'AAMI2' 'Rrhythm' 'Morph '};

    cLabListsLabels = { ...
                            cValidAnnotationsString;
                            {'Normal' 'Supraventricular' 'Ventricular' 'Fusion' 'Unclass' 'NoAAMI'};
                            {'Normal' 'Supraventricular' 'Ventricular' 'Unclass' 'NoAAMI'};
                            {'Prematuro' 'Normal' 'Escape' 'Unclass'};
                            {'Normal' 'Aberrada' 'Unclass'} };


    iLabListsTranslations = [ ...
        ...% HES labelling:    'Q' 'N' 'N' 'S' 'Q' 'V' 'S'                         
                                1:size(cTranslationMatrix,1);
                                5   1   1   2   5   3   2;
                                4   1   1   2   4   3   2;
                                4   2   2   1   4   1   1;
                                3   1   1   2   3   2   1;
                            ];

elseif( strcmpi(recording_format, 'ISHNE') )

%%%%%%%%%%%%%%%%%%
%%% ISHNE FORMAT %%%
%%%%%%%%%%%%%%%%%%%%%


    % Esta matriz es OBSOLETA. Fue reemplazada por
    % iLabListsTranslations. La deje para documentacion de las
    % clases nada mas.
    cTranslationMatrix = { ...
        'Normal beat'                                       'Normal'                1   'N'  'N' ; ...
        'Premature Ventricular Contraction'                 'Ventricular'           8   'V'  'V' ; ...
        'Supraventricular premature or ectopic beat'        'Supraventricular'      5   'S'  'S' ; ...
        'Calibration pulse'                                 'Calibration'           5   'ñ'  'C' ; ...
        'Bundel branch block beat'                          'BBB'                   9   'N'  'B' ; ...
        'Pace'                                              'Paced Beat'            4   'Q'  'P' ; ...
        'Artefact'                                          'Artefact'              9   'ñ'  'X' ; ...
        'Unknown'                                           'Unknown'               5   'Q'  'Q' ; ...
         'fv c/ NaN'                                        'fv c/ NaN'             9   'ñ'  'ñ' ; ...
                        };

    cValidAnnotations =          cTranslationMatrix(:,5);
    cAAMItranslation =           cTranslationMatrix(:,4);
    cAnnTranslation   =          cTranslationMatrix(:,2);
    cValidAnnotationsString =    cTranslationMatrix(:,1);

    cLablists = { 'ISHNE' 'AAMI' 'AAMI2' 'Rrhythm' 'Morph '};

    cLabListsLabels = { ...
                            cValidAnnotations;
                            {'Normal' 'Supraventricular' 'Ventricular' 'Fusion' 'Unclass' 'NoAAMI'};
                            {'Normal' 'Supraventricular' 'Ventricular' 'Unclass' 'NoAAMI'};
                            {'Prematuro' 'Normal' 'Escape' 'Unclass'};
                            {'Normal' 'Aberrada' 'Unclass'} };


    iLabListsTranslations = [ ...
        ...% ISHNE labelling:  'N' 'V' 'S' 'C' 'B' 'P' 'X' 'Q' 'ñ' 
                                1:size(cTranslationMatrix,1);
                                1   3   2   6   1   5   6   5   6;
                                1   3   2   5   1   3   5   4   5;
                                2   1   1   4   2   4   4   4   4;
                                1   2   1   3   2   3   3   3   3;
                            ];


elseif( strcmpi(recording_format, 'MIT') )

%%%%%%%%%%%%%%%%%
%%% MIT FORMAT %%%
%%%%%%%%%%%%%%%%%%%


    % Esta matriz es OBSOLETA. Fue reemplazada por
    % iLabListsTranslations. La deje para documentacion de las
    % clases nada mas.
    cTranslationMatrix = { ...
        'Normal beat'                                     'Normal'                                1   'N'  'N' ; ...
        'Normal beat'                                     'puntoNormal'                           1   'N'  'ñ' ; ...
        'Left bundle branch block beat'                   'L'                                     1   'N'  'L' ; ...
        'Right bundle branch block beat'                  'R'                                     1   'N'  'R' ; ...
        'Atrial premature beat'                           'APC'                                   8   'S'  'A' ; ...
        'Aberrated atrial premature beat'                 'aberrated APC'                         8   'S'  'a' ; ...
        'Nodal (junctional) premature beat'               'Nodal premature beat'                  8   'S'  'J' ; ...
        'Supraventricular premature beat'                 'S'                                     8   'S'  'S' ; ...
        'Premature ventricular contraction'               'V'                                     9   'V'  'V' ; ...
        'Fusion of ventricular and normal beat'           'F'                                     10  'F'  'F' ; ...
        'Start of ventricular flutter/fibrillation'       '['                                     11  '['  '[' ; ...
        'Ventricular flutter wave'                        'Ventricular flutter'                   11  '!'  '!' ; ...
        'End of ventricular flutter/fibrillation'         ']'                                     11  ']'  ']' ; ...
        'Atrial escape beat'                              'Atrial Escape'                         8   'S'  'e' ; ...
        'Nodal (junctional) escape beat'                  'Nodal Escape beat'                     8   'S'  'j' ; ...
        'Ventricular escape beat'                         'Ventricular Escape'                    9   'V'  'E' ; ...
        'Paced beat'                                      'Paced'                                 21  'Q'  '/' ; ...
        'Fusion of paced and normal beat'                 'fpn'                                   21  'Q'  'f' ; ...
        'Non-conducted P-wave (blocked APB)'              'Non-conducted P-wave (blocked APB)'    11  'x'  'x' ; ...
        'Non-conducted P-wave (blocked APB)'              'p'                                     11  'p'  'p' ; ...
        'Unclassifiable beat'                             'Q'                                     21  'Q'  'Q' ; ...
        'Isolated QRS-like artifact'                      'Isolated'                              11  '|'  '|' ; ...
        'Beat not classified during learning'             'Beat not classified during learning'   11  '?'  '?' ; ...
        'Rhythm change '                                  'Rhythm change'                         11  '+'  '+' ; ...
        'R-on-T premature ventricular contraction'        'R-on-T'                                9   'V'  'r' ; ...
        'ST segment change'                               'ST change'                             11  's'  's' ; ...
        'Bundle branch block beat (unspecified)'          'Bundle branch block beat'              1   'N'  'B' ; ...
        'Supraventricular escape beat (atrial or nodal)'  'Supraventricular escape beat'          8   'S'  'n' ; ...
        'fv c/ NaN'                                       'fv c/ NaN'                             11  'ñ'  'ñ' ; ...
        };

    cValidAnnotations =          cTranslationMatrix(:,5);
    cAAMItranslation =           cTranslationMatrix(:,4);
    cAnnTranslation   =          cTranslationMatrix(:,2);
    cValidAnnotationsString =    cTranslationMatrix(:,1);


    cLablists = { 'MIT-BIH' 'AAMI' 'AAMI2' 'Rrhythm' 'Morph '};

    cLabListsLabels = { ...
                            cValidAnnotationsString;
                            {'Normal' 'Supraventricular' 'Ventricular' 'Fusion' 'Unclass' 'NoAAMI'};
                            {'Normal' 'Supraventricular' 'Ventricular' 'Unclass' 'NoAAMI'};
                            {'Prematuro' 'Normal' 'Escape' 'Unclass'};
                            {'Normal' 'Aberrada' 'Unclass'} };


    iLabListsTranslations = [ ...
        ...% MITBIH labelling: 'N' 'ñ' 'L' 'R' 'A' 'a' 'J' 'S' 'V' 'F' '[' '!' ']' 'e' 'j' 'E' '/' 'f' 'x' 'p' 'Q' '|' '?' '+' 'r' 's' 'B' 'n' 'ñ'
                                1:size(cTranslationMatrix,1);
                                1   1   1   1   2   2   2   2   3   4   6   6   6   1   1   3   5   5   6   6   5   6   6   6   3   6   1   1   6;
                                1   1   1   1   2   2   2   2   3   3   5   5   5   1   1   3   4   4   5   5   4   5   5   5   3   5   1   1   5;
                                2   2   2   2   1   1   1   1   1   1   4   4   4   3   3   3   4   4   4   4   4   4   4   4   1   4   2   3   4;
                                1   1   1   1   1   2   1   1   2   2   3   3   3   1   1   2   3   3   3   3   3   3   3   3   2   3   1   1   3;
                            ];

else
    % AAMI asumption
    
    cprintf('[1,0.5,0]', 'Heartbeats annotations types in %s format no translated.\n', recording_format);
    return
                        
end

iAAMItranslation = iLabListsTranslations(2,:);

[ann_types, ~, ann_types_idx] = unique(ann.anntyp);
cant_types = length(ann_types);
[~, valid_ann_idx ann_types_valid_idx ] = intersect(cValidAnnotations, ann_types);
NotValidAnn_idx = find(iAAMItranslation(valid_ann_idx) == 6);
%descarta tambien los Unknown
% NotValidAnn_idx = find(iAAMItranslation(valid_ann_idx) >= 5);
NotValidAnn_idx = ann_types_valid_idx(NotValidAnn_idx);

%Termino de averiguar las anotaciones que no sean latidos.
aux = setdiff(1:cant_types, ann_types_valid_idx);
NotValidAnn_idx = [(NotValidAnn_idx(:))' aux];
ValidAnn_idx = setdiff(1:length(ann_types), NotValidAnn_idx);

[~, valid_ann_idx ] = intersect(cValidAnnotations, ann_types(ValidAnn_idx));

ann_types_translated = zeros(cant_types,1);
ann_types_translated(ValidAnn_idx) = valid_ann_idx;

%Resguardar� un vector con el c�digo de anotacion para cada latido.
ann_types_idx = ann_types_translated(ann_types_idx)';

bValidBeats = colvec(ann_types_idx ~= 0);

for field = rowvec(fieldnames(ann))
    if( isempty(ann.(field{1})) )
        ann = rmfield(ann, field{1} );
    else
        ann.(field{1}) = ann.(field{1})(bValidBeats);
    end
end

%Convert to the proper labeling
ann.anntyp = char(cAAMItranslation(colvec(iLabListsTranslations( LabConversion_idx, ann_types_idx(bValidBeats)))));
