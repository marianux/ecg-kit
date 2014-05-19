function bFileFound = exist_distributed_file(file_name, retries)

bFileFound = false;
while( ~bFileFound && retries >= 0)
    if( exist(file_name, 'file') )
        bFileFound = true;
    else            
       retries = retries - 1;
       pause(10)
    end
end
