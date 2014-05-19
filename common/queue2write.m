return

if( ~bHaveUserInterface )
    jj = Loops2io;
    while(jj>0)
        loc_files = dir([tmp_path '*.wloc']);
        lloc_files = length(loc_files);
        if( lloc_files < MaxNodesWriting  )
            jj = 0;
        else
            jj = jj - 1;
            pause((5+2*rand(1))); 
        end
    end
    
    %tomo posecion del medio.
    fid_loc = fopen( [tmp_path rec_filename '_' num2str(this_pid) '.wloc'] , 'w');
    fclose(fid_loc);
    
end
