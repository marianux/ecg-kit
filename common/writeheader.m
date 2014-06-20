function writeheader(header_path, header)
%
% WRITHEAD function writes header file for signal data struct in directory header_path
%	Input parameters.
%	   header_path: directory of work
%	   header: struct of header to write

if(header_path(end) ~= filesep )
    header_path = [header_path filesep];
end

% Opening ASCII header file
if ~isempty(header.recname)   
    fid=fopen([header_path  header.recname '.hea'], 'w');
    header.fname = repmat([ header.recname '.dat'], header.nsig,1);
else
    fid = fopen([header_path  header.fname(1,:) '.hea'], 'w');
end

if (fid<0)
    error(['Can not write ' [header_path  header.recname '.hea'] '\n']);
end

try

    % writing first line of record_name, # signals and so on
    fprintf(fid,'%s %d %d %d\n',header.recname,header.nsig);
    if isfield(header,'freq')
        fprintf(fid,'%d ',header.freq);
    end
    if isfield(header,'nsamp')
        fprintf(fid,'%d ',header.nsamp);
    end
    if isfield(header,'btime')
        fprintf(fid,'%s ',header.btime);
    end
    if isfield(header,'bdate');
        fprintf(fid,'%s ',header.bdate);
    end

    fmt = 16;
    % writing the rest of lines
    for ii=1:header.nsig

        if isfield(header,'spf')||isfield(header,'skew')||isfield(header,'offset')
            if isfield(header,'spf')
                fprintf(fid,'\n%s %d',header.fname(ii,:),fmt);
                fprintf(fid,'x%d ',header.spf(ii));
            end
            if isfield(header,'skew')
                fprintf(fid,'\n%s %d',header.fname(ii,:),fmt);
                fprintf(fid,':%d ',header.skew(ii));
            end
            if isfield(header,'offset')
                fprintf(fid,'\n%s %d',header.fname(ii,:),fmt);
                fprintf(fid,'+%d ',header.offset(ii));
            end
        else
            fprintf(fid,'\n%s %d ',header.fname(ii,:),fmt);
        end

        if isfield(header,'gain')
            fprintf(fid,'%f',header.gain(ii) );
        else
            fprintf(fid,'%f',200 );
        end
        
        if isfield(header,'baseline')
            if isfield(header,'units')
                fprintf(fid,'(%d)',header.baseline(ii));
            else
                fprintf(fid,'(%d) ',header.baseline(ii));
            end
        elseif (isfield(header,'units') )
            fprintf(fid,'/%s ',header.units(ii,:));
        else
            fprintf(fid,' ' );
        end

        if isfield(header,'adcres')
            fprintf(fid,'%d ',header.adcres(ii));
        else
            fprintf(fid,'%d ', 12);
        end
        
        if isfield(header,'adczero')
            fprintf(fid,'%d ',header.adczero(ii));
        else
            fprintf(fid,'%d ', 0);
        end
        
        if isfield(header,'initval')
            fprintf(fid,'%d ',header.adczero(ii));
        else
            fprintf(fid,'%d ', 0);
        end
        
        if isfield(header,'cksum')
            fprintf(fid,'%d ',header.cksum(ii));
        else
            fprintf(fid,'%d ', 0);
        end
        
        if isfield(header,'bsize')
            fprintf(fid,'%d ',header.bsize(ii));
        else
            fprintf(fid,'%d ', 0);
        end
        
        if isfield(header,'desc')
            fprintf(fid,'%s',header.desc(ii,:));
        end    
        
    end

    fprintf(fid,'\n');
    
    fclose(fid);

catch MEE
    fclose(fid);
    rethrow(MEE);
end
