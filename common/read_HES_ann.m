%% Reads ECG annotations in HES format
% Reads ECG annotations in HES (Biosigna) format. Implements the documentation
% available in the help of the application provided with the database (not
% available with the ECHkit). 
% 
% Arguments:
%   + filename: recording to be read.
% 
% Output:
%   + ann: annotations for the ECG recordings.
% 
% Limits:
% This routine is limited to read blocks smaller than 200 Mbytes for
% performance reasons. You can disable this limit by doing:
% MaxIOread = Inf; %megabytes
% 
% See also read_HES_format, read_HES_header, read_ECG, ECGwrapper
% 
% Author: Mariano Llamedo Soria
% <matlab:web('mailto:llamedom@electron.frba.utn.edu.ar','-browser') (email)> 
% Version: 0.1 beta
% Birthdate: 17/12/2010
% Last update: 19/11/2014
% Copyright 2008-2014
% 
function ann = read_HES_ann(filename)

% La traduccion se termina haciendo en LoadDatabaseCustomization.m
% beat_type_idx = [ 0 1 2 4 5 60 70 ];
% beat_type_desc = { ...
%                     'NO_TYPING' ...
%                     'DOMINANT_TYPE_0' ...
%                     'DOMINANT_TYPE_1' ...
%                     'ABERRANT' ...
%                     'ARTEFAKT' ...
%                     'VES' ...
%                     'SVES' ...
%                 };
% beat_type_AAMI_translation = [ 'Q';'N';'N';'S';'Q';'V';'S' ];

fid = fopen(filename);

if( fid > 0 )

    bDesktop = usejava('desktop') & ispc();
    
    try 

        bDataNotFound = true;
        
        iLineCount = 0;
        
        while( ~feof(fid) && bDataNotFound && iLineCount < 20 )
        
            strLine = fgetl(fid);
            iLineCount = iLineCount + 1;
            
            [dummy tokens]= regexp(strLine, '.*(\d+):(\d+):(\d+):(\d+).*', 'match', 'tokens');

            if( ~isempty(tokens) )
                bDataNotFound = false;
            end
                
        end
        
        if( ~bDataNotFound )
            
            frewind(fid);            
            % Se encontro la tabla de anotaciones -> leemos todo lo demas
            lines = textscan(fid, '%s', 'delimiter', '\n');

            fclose(fid);
        end
       
    catch any_err
        fclose(fid)
        rethrow(any_err)
    end

    try
        
        if( ~bDataNotFound )    
            iBeatCount = 1;
            iArtifactCount = 1;

            lines = lines{1};
            lines = lines(iLineCount+1:2:end);
            llines = length(lines);

            ann.time = nan(llines,1);
            ann.anntyp = nan(llines,1);
            
            if(bDesktop)            
                % Activate the progress_struct bar.
                pb = progress_bar(filename);
            end
            ii = 0;

            pb.Loops2Do = length(lines);
            
            for cell_line = lines'

                if(bDesktop)
                    %start of the progress_struct loop 0%
                    pb.start_loop();
                end
                
                strLine = cell_line{1};
                
                [dummy tokens]= regexp(strLine, '\s*(\d+)\s+(\d+):(\d+):(\d+):(\d+)\s+(\d+)\s+.*', 'match', 'tokens');

                if( ~isempty(tokens) )
                    
                    tokens = tokens{1};

                    BeatNum = ceil(str2num(tokens{1})/2);

                    if(iBeatCount ~= BeatNum)
                       error('Problema de sincronizacion en la lectura, revisar!!') 
                    end

                    BeatTime = str2num(tokens{2})*60*60 + str2num(tokens{3})*60 + str2num(tokens{4}) + str2num(tokens{5}) * 1e-3;
                    ann.time(BeatNum) = BeatTime;
                    
                    BeatType = str2num(tokens{6});
                    
                    ann.anntyp(BeatNum) = BeatType;
                    
                    iBeatCount = iBeatCount + 1;

                end

                ii = ii + 1;
                
                if(bDesktop)
                    pb.checkpoint('parsing annotations ...');
                    
                    pb.end_loop();
                end

                
            end

            if(bDesktop)            
                % destroy the progress bar.
                pb.delete;
            end

            ann.time = ann.time(1:iBeatCount-1);
            ann.anntyp = ann.anntyp(1:iBeatCount-1);
            
        else
            error(['No se encontró el patrón de sincronización en ' filename])
        end

    catch any_err
        
        if(bDesktop)
            if( ~isempty(wb_h))
                % destroy the progress bar.
                pb.delete;
            end
        end
        
        rethrow(any_err);
        
    end

    
else
    error('No se pudo abrir el archivo..')
end
