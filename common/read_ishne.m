
function [ ECG heasig ann last_sample ] = read_ishne(filename, start_sample, end_sample)

% Reads ECG recordings in ISHNE format from THEW databases. Implements the documentation available in:
% 
% http://thew-project.org/THEWFileFormat.html
% 
% Arguments:
%   + filename: recording to be read.
%   + start_sample: (opt) start sample to read. Default 1.
%   + end_sample: (opt) end sample to read. Default min(All recording, ECG block of 200 Mbytes)
% 
% Output:
%   + ECG: the ECG block
%   + heasig: header with the ECG properties. 
%   + ann: annotations for the ECG recordings.
% 
% Limits:
% This routine is limited to read blocks smaller than 200 Mbytes for
% performance reasons. You can disable this limit by doing:
% MaxIOread = Inf; %megabytes
% 
% Author: Mariano Llamedo Soria (llamedom at {unizar.es;electron.frba.utn.edu.ar}
% Birthdate: 21/7/2010
% Last update: 20/02/2013
% 

%No leer bloques mas grandes de 200 megabytes
MaxIOread = 200; %megabytes

if( nargin < 2 || isempty( start_sample ) )
    start_sample = 1;
else
    start_sample = max(1,start_sample);
end

ann = [];
heasig = [];
ECG = [];

if( nargout > 1 )
    bHeaderRequired = true;
else
    bHeaderRequired = false;
end

if( nargout > 2 )
    bAnnRequired = true;
    ann_filename = [filename(1:end-4) '.ann'];
else
    bAnnRequired = false;
    ann_filename = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%
% primero leo el header.
    
% if(bHeaderRequired)
    heasig = read_ishne_header(filename);
% end


%%%%%%%%%%%%%%%%%%%%%%
% Ahora leo la señal.

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

        % Corroboro la integridad header vs contenido
        fseek(fid, 0, 'eof');
        bytes_totales = ftell(fid);
        muestras_guardadas = (bytes_totales - (522+var_length_size)) / 2 / heasig.nsig;
        
        if( abs(muestras_guardadas - heasig.nsamp) > 584 )
            heasig.nsamp = muestras_guardadas;
            warning('read_ishne:CorruptFile', 'File corrupted. Check data/header integrity.')
        end
        
        bContinue = true;
        if( nargin < 3 || isempty( end_sample ) )
            %Intento la lectura total por defecto
            samples2read = heasig.nsamp - (start_sample-1);
        else
            samples2read = min(heasig.nsamp, end_sample) - (start_sample-1);
        end

        if( (samples2read*heasig.nsig*2) > (MaxIOread * 1024^2) )
            samples2read = (MaxIOread * 1024^2) / heasig.nsig / 2;
            warning(['No es recomendable leer mas de ' num2str(MaxIOread) ' Mb. Realice varias lecturas.'])
        end

        if(samples2read <= 0)
            error('read_ishne:Impossible2read', ['Imposible leer ' num2str(samples2read) ' muestras'])
        end

        bWarn_menos_muestras = false;

%         if( ispc )
%             %Desactivo esto para que no me moleste ni al debuguear.
%             dbclear if caught error	
%         end

        while( bContinue )

            try

                fseek(fid, 522+var_length_size+((start_sample-1)*heasig.nsig)*2, 'bof');
                ECG = fread(fid, [heasig.nsig samples2read], '*int16')';
                bContinue = false;

            catch ME
                bWarn_menos_muestras = true;
                samples2read = round(samples2read / 2);
            end

        end

%         if( ispc )
%             dbstop if caught error	
%         end

        fclose(fid);
    
    catch ME2
        
        fclose(fid);
        
        rethrow(ME2)
        
    end
    
    samples_actually_read = size(ECG,1);
    
    last_sample = samples_actually_read + start_sample - 1;
    
    if( samples_actually_read < samples2read )
        warning('Fin del archivo alcanzado antes de los esperado. Corroborar si es correcto el header, o el registro está truncado.');
    end
    
    if( bWarn_menos_muestras ) 
        warning(['Límite de memoria del sistema excedido. Se leyó de ' num2str(start_sample) ':' num2str(last_sample)])
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ahora leo las anotaciones si las hay.

if( bAnnRequired )
    ann = read_ishne_ann(ann_filename);
end
