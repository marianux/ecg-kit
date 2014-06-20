function wfdb_prefix = WFDB_command_prefix(recording_path)

    % string to change the working directory.
    if( ispc() )
        % in windows, it is important to change the
        % drive unit before the directory change
        % command (CD)
        command_sep_str = '&';
        wfdb_prefix = recording_path;
        wfdb_prefix = [wfdb_prefix(1:2) command_sep_str ];
    else
        command_sep_str = ';';
        wfdb_prefix = [];
    end
    wfdb_prefix = [wfdb_prefix 'cd ' recording_path command_sep_str ];
