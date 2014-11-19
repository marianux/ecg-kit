%% Reads ECG header from ISHNE format
% Reads ECG recordings in ISHNE format from THEW databases. Implements the
% documentation available in: 
% 
% http://thew-project.org/THEWFileFormat.html
% 
% Arguments:
%   + filename: recording to extract header block.
% 
% Output:
%   + heasig: header with the ECG properties. 
% 
% 
% See also read_ishne_ann, read_ishne_format, read_ECG, ECGwrapper
% 
% Author: Mariano Llamedo Soria
% <matlab:web('mailto:llamedom@electron.frba.utn.edu.ar','-browser') (email)> 
% Version: 0.1 beta
% Birthdate: 21/7/2010
% Last update: 19/11/2014
% Copyright 2008-2014
% 
function heasig = read_ishne_header(filename)

lead_transformation =   { ...
                        'Unknown'; ...
                        'bipolar'; ...
                        'X'; ...
                        'Y'; ...
                        'Z'; ...
                        'I'; ...
                        'II'; ...
                        'III'; ...
                        'aVL'; ... ¿ Esto no se si está bien o al reves?
                        'aVR'; ... ¿ Esto no se si está bien o al reves?
                        'aVF'; ...
                        'V1'; ...
                        'V2'; ...
                        'V3'; ...
                        'V4'; ...
                        'V5'; ...
                        'V6'; ...
                        'ES'; ...
                        'AS'; ...
                        'AI'; ...
                        };

heasig = [];

fid = fopen(filename);

if( fid > 0 )

    try 
        magic_num = fread(fid, 8,'*char');

        if( feof(fid) || ~strcmpi(rowvec(magic_num(1:5)), 'ISHNE') )
            fclose(fid);
            return
        end

        %salteo el crc
        fseek(fid, 2, 'cof');

        var_length_size = fread(fid, 1,'int32');

        heasig.nsamp = fread(fid, 1,'int32');

        %Paso a la fecha
        fseek(fid, 120, 'cof');
        date_of_recording = fread(fid, 3,'int16');
        date_of_file_creation = fread(fid, 3,'int16');
        time_of_recording = fread(fid, 3,'int16');
        heasig.btime = sprintf('%0d:%0d:%0d',time_of_recording(1),time_of_recording(2),time_of_recording(3));
        heasig.bdate = sprintf('%0d/%0d/%4d',date_of_recording(1),date_of_recording(2),date_of_recording(3));
        heasig.nsig = fread(fid, 1,'int16');

%         al parecer esto ya no es más así
%         %nsamp era la cantidad total de muestras en todos los canales.
%         heasig.nsamp = heasig.nsamp / heasig.nsig;

        lead_spec = fread(fid, 12,'int16');

        heasig.desc = char(lead_transformation(lead_spec(1:heasig.nsig)+1));

        lead_quality = fread(fid, 12,'int16');

        lead_amp_res = fread(fid, 12,'int16');

        heasig.gain = 1./(lead_amp_res(1:heasig.nsig));

        pacemaker_code = fread(fid, 1,'int16');

        recorder_type = char(fread(fid, 40,'uchar'));

        heasig.freq = fread(fid, 1,'int16');

        heasig.adcres = repmat(16,heasig.nsig,1);
        heasig.adczero = zeros(heasig.nsig,1);
        heasig.units = repmat('nV',heasig.nsig,1);
    
        [~, recname] = fileparts(filename);
        
        heasig.recname = recname;
        
        % Corroboro la integridad header vs contenido
        fseek(fid, 0, 'eof');
        bytes_totales = ftell(fid);
        muestras_guardadas = round((bytes_totales - (522+var_length_size)) / 2 / heasig.nsig);
        
        if( abs(muestras_guardadas - heasig.nsamp) > 584 )
            str_aux = sprintf('Header of %s corrupted, expected %d samples and found %d. Reading only the available samples. Check data/header integrity.', heasig.recname, heasig.nsamp, muestras_guardadas);
            bUseDesktop = usejava('desktop');
            if(bUseDesktop)
                fprintf(1, ' ')
                cprintf('[1,0.5,0]', str_aux);
                fprintf(1, '\n')
            else
                % esto lo quito porque ensucia los logs.
%                 fprintf(2, '%s\n', str_aux);
            end
            heasig.nsamp = muestras_guardadas;
        end
        
        fclose(fid);
        
    catch ME
        fclose(fid);
        
        rethrow(ME)
    end
    
end
