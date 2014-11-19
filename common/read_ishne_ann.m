%% Reads ECG annotations from ISHNE format
% Reads ECG recordings in ISHNE format from THEW databases. Implements the
% documentation available in: 
% 
% http://thew-project.org/THEWFileFormat.html
% 
% Arguments:
%   + ann_filename: annotation file to read annotations.
% 
% Output:
%   + ann: structure of annotations. 
% 
% See also read_ishne_ann, read_ishne_header, read_ECG, ECGwrapper
% 
% Author: Mariano Llamedo Soria
% <matlab:web('mailto:llamedom@electron.frba.utn.edu.ar','-browser') (email)> 
% Version: 0.1 beta
% Birthdate: 21/7/2010
% Last update: 19/11/2014
% Copyright 2008-2014
% 
function ann = read_ishne_ann( ann_filename )

ann = [];
fid = fopen(ann_filename);

if( fid > 0 )

    try
        magic_num = fread(fid, 8,'*char');

        if( feof(fid) || ~strcmpi(rowvec(magic_num(1:3)), 'ann') )
            fclose(fid);
            return
        end

        %salteo el crc
        fseek(fid, 2, 'cof');

        var_length_size = fread(fid, 1,'int32');

        %voy al comienzo de las anotaciones
        if( feof(fid) )
            fclose(fid);
            return
        end
        fseek(fid, 522+var_length_size, 'bof');

        init_sample = fread(fid, 1,'uint32');

        fseek(fid, 0, 'eof');
        filesize = ftell(fid);

        ann_size = (filesize - (526+var_length_size))/4;

        fseek(fid, 526+var_length_size, 'bof');

        lab1 = char(fread(fid, ann_size, 'uint8', 3));

        fseek(fid, 527+var_length_size, 'bof');

        lab2 = char(fread(fid, ann_size, 'uint8', 3));

        fseek(fid, 528+var_length_size, 'bof');

        toc =  fread(fid, ann_size, 'uint16', 2);

        fclose(fid);
    
    catch ME
        fclose(fid);
        
        rethrow(ME)
    end

    bTimeouts = lab1 == '!';

    toc = init_sample + cumsum(toc);

    %Quito las anotaciones de timeouts.
    toc(bTimeouts) = [];
    lab1(bTimeouts) = [];
    lab2(bTimeouts) = [];

    ann.time = toc;
    ann.anntyp = lab1;
    ann.subtyp = lab2;

end
