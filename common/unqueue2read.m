return

if( ~bHaveUserInterface )
    loc_file = dir([tmp_path rec_filename '_' num2str(this_pid) '.rloc']);
    if( ~isempty(loc_file) )
        %tomo posecion del medio.
        delete([tmp_path rec_filename '_' num2str(this_pid) '.rloc']);
    end
end
