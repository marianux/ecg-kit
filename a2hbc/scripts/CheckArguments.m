
% op_mode parsing
if( isnumeric(op_mode) )
    op_mode = cKnownModesOfOperation{op_mode};
end

% class filtering parsing
class_filter = char(intersect( cellstr(class_filter), cellstr(typical_lablist)));
if( isempty(class_filter) )
    strAux = [ repmat(' + ', size(typical_lablist,1), 1) char(typical_lablist) repmat('\n', size(typical_lablist,1), 1 ) ];
    error( 'a2hbc:ArgCheck:InvalidClassFilter', ['Invalid class-filter. Please provide one of these classes:\n' rowvec(strAux')] );
end

%ECG parsing
if( ~isempty(ECG_total) )
    % ECG already read

    if( isempty(ECG_header))
       error( 'a2hbc:ArgCheck:InvalidHeader', 'Please provide the ECG header.\n\n' );
    else
        if( ~isfield(ECG_header, cHeaderFieldNamesRequired ) )
            strAux = [ repmat(' + ', length(cHeaderFieldNamesRequired), 1) char(cHeaderFieldNamesRequired) repmat('\n', length(cHeaderFieldNamesRequired), 1 ) ];
            error( 'a2hbc:ArgCheck:InvalidHeader', ['Please provide the following fields in the header struct:\n ' rowvec(strAux') ] );
        end
    end
    
    heasig = ECG_header;
    
    if( isempty(ECG_annotations))
       error( 'a2hbc:ArgCheck:InvalidHeader', 'Please provide the ECG annotations.\n\n' );
    else
        if( ~isfield(ECG_annotations, cAnnotationsFieldNamesRequired ) )
            strAux = [ repmat(' + ', length(cAnnotationsFieldNamesRequired), 1) char(cAnnotationsFieldNamesRequired) repmat('\n', length(cAnnotationsFieldNamesRequired), 1 ) ];
            error( 'a2hbc:ArgCheck:InvalidAnnotations', ['Please provide the following fields in the annotations struct:\n ' rowvec(strAux') ] );
        end
    end
    
    ann = ECG_annotations;
    
    [QRS_locations, true_labels ] = Annotation_process(ann, recording_format, class_labeling);
    
elseif( ~isempty(recording_name) )
    
    [~, rec_filename] = fileparts(recording_name);
    
    % ECG to be read, create a wrapper object, later assign tasks to it.
    ECGw = ECGwrapper( ... 
                        'recording_name', recording_name, ...
                        'recording_format', recording_format, ...
                        'this_pid', this_pid, ...
                        'cant_pids', cant_pids ...
                        );
    
    try
        ECGw.CheckECGrecording();
    catch ME
        thisME = MException( 'a2hbc:ArgCheck:InvalidECGarg', 'Please provide an ECG recording as described in the documentation, help(''a2hbc'') maybe could help you.\n' );
        ME = addCause(ME, thisME );
        rethrow(ME);
    end
    
else
%     strAux = help('a2hbc'); %#ok<MCHLP>
    error( 'a2hbc:ArgCheck:InvalidECGarg', 'Please provide an ECG recording as described in the documentation, help(''a2hbc'') maybe could help you.\n' );
end

bLabelingChanged = false;

lablist_idx = find(strcmpi(class_labeling, cKnownLabelings));

%Update the automatic classifier.
global_classifier = global_classifiers{lablist_idx};
                
if( isempty(tmp_path) )
    tmp_path = tempdir;
end

%check path integrity.
if(~exist(tmp_path, 'dir'))
    %try to create it
    if( ~mkdir(tmp_path) )
        error('a2hbc:ArgCheck:InvalidPath', 'Invalid tmp_path. Please provide a valid path.\n' );
    end
end

bArgCheckPassed = true;
           
