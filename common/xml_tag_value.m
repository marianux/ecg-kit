%% (Internal) Checks (tagName, tagValue) in an XML node "element"
% 
% 
% See also read_hl7a_format
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 09/05/2016
% Birthdate  : 09/05/2016
% Copyright 2008-2016
% 
function [bRetVal, tagActualValue] = xml_tag_value(element, tagName, tagValue)

    bRetVal = false;
    tagActualValue = [];
    
    if( nargin < 3 )
        tagValue = '';
    end

    thisCode_att = element.getAttributes;        

    for ii = 0:(thisCode_att.getLength-1)

        this_att = thisCode_att.item(ii);

        aux_val = this_att.getName;

        if( strcmpi(aux_val, tagName) )
            tagActualValue = char(this_att.getValue);
            bRetVal = strcmpi(tagActualValue, tagValue);
            return;
        end
        
    end

    