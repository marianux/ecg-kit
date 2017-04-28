%% (Internal) Write an ECG header in MIT format
% 
% This function receives a "header" struct with ECG signal propeties and
% writes [header.name '.hea'] file into "header_path".
% 
%       writeheader(header_path, header)
% 
% Arguments:
% 
%	   header_path: directory of work
% 
%	   header: struct of header to write
% 
% Example
% 
%         MIT_filename  = 'your_filename';
%         MIT_path = ['.' filesep];
%         ECG_header.recname = MIT_filename;
%         ECG_header.nsig = 8;
%         ECG_header.nsamp = 9978;
%         ECG_header.freq = 1000; %Hz
%         ECG_header.desc = char({'I','II','V1','V2','V3','V4','V5','V6'});
%         ECG_header.adczero = zeros(ECG_header.nsig,1);
%         ECG_header.gain = ones(ECG_header.nsig,1);
%         ECG_header.units = repmat('uV',ECG_header.nsig,1);
%         ECG_header.btime = '00:00:00';
%         ECG_header.bdate = '01/01/2000';
% 
%         writeheader(MIT_path, ECG_header);
% 
% See also writeannot, read_ECG, ECGwrapper
% 
% Author: Salvador Olmos, Juan Pablo Mart�nez and/or Juan Bolea
% adapted to ecg-kit by: Mariano Llamedo Soria
% <matlab:web('mailto:llamedom@electron.frba.utn.edu.ar','-browser') (email)> 
% Version: 0.1 beta
% Birthdate: 5/01/2014
% Last update: 19/11/2014
% Copyright 2008-2015
% 
function writeheader(header_path, header)
%

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
    error('Can not write %s.hea\n', header_path, header.recname );
end

try

    % writing first line of record_name, # signals and so on
    fprintf(fid,'%s %d %d %d\n',header.recname, header.nsig);
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

        if isfield(header,'spf') || isfield(header,'skew') || isfield(header,'offset')
            fprintf(fid,'\n%s %d',header.fname(ii,:),fmt);
            if isfield(header,'spf') && ~isnan(header.spf(ii))
                fprintf(fid,'x%d',header.spf(ii));
            end
            if isfield(header,'skew') && ~isnan(header.skew(ii))
                fprintf(fid,':%d',header.skew(ii));
            end 
            if isfield(header,'offset') && ~isnan(header.offset(ii))
                fprintf(fid,'+%d',header.offset(ii));
            end
            fprintf(fid,' ');
        else
            fprintf(fid,'\n%s %d ',header.fname(ii,:),fmt);
        end

        if isfield(header,'gain')
            fprintf(fid,'%f',header.gain(ii) );
        else
            fprintf(fid,'%f',200 );
        end
        
        if(isfield(header,'baseline') && ~isnan(header.baseline(ii)) ) 
            fprintf(fid,'(%d)',header.baseline(ii));
        end
        
        if (isfield(header,'units') )
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
            fprintf(fid,'%d ',header.initval(ii));
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
