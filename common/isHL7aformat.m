%% (Internal) Check if a recording is in HL7a format.
%   
%   bRetval = isHL7aformat(filename)
% 
% Arguments:
% 
%      + filename: the recording
% 
% Output:
% 
%      + bRetval: Boolean if it is of this format.
% 
% Example:
% 
% See also ECGformat, read_ECG, isAHAformat
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Birthdate: 04/05/2016
% Last update: 04/05/2016
% Copyright 2008-2016
% 
function bRetval = isHL7aformat(filename)

bRetval = false;

try
   xDoc = xmlread(filename);
catch
%     error('isHL7aformat:ReadError', 'Failed to read XML file %s.\n', filename);
    return
end

AnnECGtag = xDoc.getElementsByTagName('AnnotatedECG');

if( AnnECGtag.getLength > 0 )

    Codetag = xDoc.getElementsByTagName('code');
    
    thisVal = Codetag.item(0);

    thisVal_att = thisVal.getAttributes;

    for ii = 0:(thisVal_att.getLength-1)

        this_att = thisVal_att.item(ii);

        aux_val = char(this_att.getName);

        if( strcmpi(aux_val, 'code') )

            if( strcmpi( char(this_att.getValue), '93000' ) )
                % HL7a format fingerprint
                bRetval = true;
                return;
            end
            
        end

    end

end


