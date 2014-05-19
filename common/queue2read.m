
return

if( ~bHaveUserInterface )
    jj = Loops2io;
    while(jj>0)
        loc_files = dir([tmp_path '*.rloc']);
        lloc_files = length(loc_files);
        if( lloc_files < MaxNodesReading  )
            %tomo posecion del medio.
            jj = 0;
        else
            jj = jj - 1;
            pause((5+2*rand(1))); 
        end
    end
    
    fid_loc = fopen( [tmp_path rec_filename '_' num2str(this_pid) '.rloc'] , 'w');
    fclose(fid_loc);
    
end
