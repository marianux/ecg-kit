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
wb_h = [];

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
                wb_h = waitbar(0,'Parsing annotations...');
            end
            ii = 0;

            for cell_line = lines'

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
                    waitbar(ii/llines, wb_h);
                end

            end

            if(bDesktop)            
                close(wb_h)
            end

            ann.time = ann.time(1:iBeatCount-1);
            ann.anntyp = ann.anntyp(1:iBeatCount-1);
            
        else
            error(['No se encontró el patrón de sincronización en ' filename])
        end

    catch any_err
        
        if(bDesktop)
            if( ~isempty(wb_h))
                close(wb_h);
            end
        end
        
        rethrow(any_err);
        
    end

    
else
    error('No se pudo abrir el archivo..')
end
