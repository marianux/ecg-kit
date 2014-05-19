return

if( ~bHaveUserInterface )
    loc_file = dir([tmp_path rec_filename '_' num2str(this_pid) '.wloc']);
    if( ~isempty(loc_file) )
        %tomo posecion del medio.
        delete([tmp_path rec_filename '_' num2str(this_pid) '.wloc']);
    end
end
