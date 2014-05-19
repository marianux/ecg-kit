function heasig = read_ishne_header(filename)

% Reads ECG headers in ISHNE format from THEW databases. Implements the documentation available in:
% 
% http://thew-project.org/THEWFileFormat.html
% 
% Arguments:
%   + filename: recording to extract header block.
% 
% Output:
%   + heasig: header with the ECG properties. 
% 
% Limits:
% This routine identifies ECG leads as described in the documentation
% available in the THEW web page. Whereas, the author is not sure after
% reading this documentation that the leads aVL and aVR are identified
% correctly (See lead_transformation). Any feedback is very welcome.
% 
% Author: Mariano Llamedo Soria (llamedom at {unizar.es;electron.frba.utn.edu.ar}
% Birthdate: 21/7/2010
% Last update: 18/2/2013
% 


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
        date_of_creation = fread(fid, 3,'int16');
        time_of_recording = fread(fid, 3,'int16');
        heasig.date = datenum(date_of_recording(3),date_of_recording(2),date_of_recording(1),time_of_recording(1),time_of_recording(2),time_of_recording(3));

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
        muestras_guardadas = (bytes_totales - (522+var_length_size)) / 2 / heasig.nsig;
        
        if( abs(muestras_guardadas - heasig.nsamp) > 584 )
            heasig.nsamp = muestras_guardadas;
            warning('read_ishne:CorruptFile', 'File corrupted. Check data/header integrity.')
        end
        
        fclose(fid);
        
    catch ME
        fclose(fid);
        
        rethrow(ME)
    end
    
end
