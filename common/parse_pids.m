function [this_pid, cant_pids] = parse_pids( aux_str )

if( ischar(aux_str) )
     delim = ' /-';
    [this_pid,aux_str] = strtok(aux_str,delim);
    cant_pids = strtok(aux_str, delim);
    this_pid = str2double(this_pid);
    cant_pids = str2double(cant_pids);
else
    if( length(aux_str) > 1 )
        this_pid = aux_str(1);
        cant_pids = aux_str(2);
    else
        this_pid = aux_str(1);
        cant_pids = aux_str(1);
    end
end
