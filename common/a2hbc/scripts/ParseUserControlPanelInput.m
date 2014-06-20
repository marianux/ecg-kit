
if( ishandle(UCP_struct.fig_hdl) )
    %classifier arguments
    aux_str = get(UCP_struct.recording_name, 'String');
     
    %Check if reprocessing is needed
    if( ~strcmpi(recording_name, aux_str) )
        bRecordingChanged = true;
    end
    recording_name = aux_str;
    
    recording_format = cKnownFormats{get(UCP_struct.recording_fmt, 'Value')};
    tmp_path = get(UCP_struct.tmp_path, 'String');
    op_mode = cKnownModesOfOperation{get(UCP_struct.op_mode, 'Value')};

    %Labeling options
    aux_str = cellstr(get(UCP_struct.labelings, 'String'));
    aux_val = get(UCP_struct.labelings, 'value');
     
    %Check if reprocessing is needed
    if( ~strcmpi(class_labeling, aux_str{aux_val}) )
        bLabelingChanged = true;
    end
    class_labeling = aux_str{aux_val};
    
    aux_str = cellstr(get(UCP_struct.class_filter, 'string'));
    aux_val = get(UCP_struct.class_filter, 'value');
    
    class_filter = char(aux_str(aux_val));
    
    %clustering options
    CantClusters = get(UCP_struct.CantClusters, 'Value');
    iter_times = get(UCP_struct.iter_times, 'Value');
    cluster_presence = get(UCP_struct.cluster_presence, 'Value');

    %User interface options
    CantSamplesClose = str2double(get(UCP_struct.CantSamplesClose, 'String'));
    CantSamplesFar = str2double(get(UCP_struct.CantSamplesFar, 'String'));
end
