%% (Internal) Display as formated strings an XML node array
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
function strAux = disp_xml_nodes(thisNode)

    strAux = [];
    
    for jj = 0:(thisNode.getLength-1)
    
        element = thisNode.item(jj);

        strAux = sprintf( '%s <%s' , strAux, char(element.getTagName) );
        
        thisCode_att = element.getAttributes;        
        
        for ii = 0:(thisCode_att.getLength-1)

            this_att = thisCode_att.item(ii);

            strAux = sprintf( '%s %s ' , strAux, char(this_att.toString) );

        end
        
        strAux = sprintf( '%s>\n' , strAux );
        
    end
        