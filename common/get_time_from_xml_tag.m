%% (Internal) Converts time from tag to seconds in an XML node "element"
% 
% 
% See also read_hl7a_format
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 10/05/2016
% Birthdate  : 10/05/2016
% Copyright 2008-2016
% 
function time_in_seconds = get_time_from_xml_tag(element, bTimeRelative)

    if( nargin < 2 || ~islogical(bTimeRelative) )
        bTimeRelative = true;
    end

    time_in_seconds = [];
    aux_time = 1;
    bReadOk = 0;

    if( isempty(element) ) 
        return
    end
    
    thisCode_att = element.getAttributes;        

    for ii = 0:(thisCode_att.getLength-1)

        this_att = thisCode_att.item(ii);

        aux_val = this_att.getName;

        if( bTimeRelative )
            
            if( strcmpi(aux_val, 'value') )
                aux_time = aux_time * str2double(char(this_att.getValue));
                bReadOk = bReadOk + 1;
            elseif( strcmpi(aux_val, 'unit') )
                aux_val = char(this_att.getValue);
                bReadOk = bReadOk + 1;
                switch(aux_val)
                    case 's'
                        aux_time = aux_time;
                    case 'ms'
                        aux_time = aux_time * 1e-3;
                    case 'us'
                        aux_time = aux_time * 1e-6;
                    otherwise
                        error('get_time_from_xml_tag:ParseError', 'Parse error, check unit = %s\n', aux_val);
                end

            end
            
        else
            % time absolute
            if( strcmpi(aux_val, 'value') )
                aux_val = char(this_att.getValue);
                
                if( length(aux_val) == 18 )
                    aux_time = datenum(aux_val,'yyyymmddHHMMSS.FFF' );
                    bReadOk = 2;
                elseif( length(aux_val) == 14 )
                    aux_time = datenum(aux_val,'yyyymmddHHMMSS' );
                    bReadOk = 2;
                end
                
            end
            
        end
        if( bReadOk == 2 )
            time_in_seconds = aux_time;
            return
        end
        
    end
    
    