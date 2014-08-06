function [strRetVal, iHours, iMins, iSeconds, iMilliSeconds ] = Seconds2HMS(data, prec)

if( nargin < 2 )
    % decimals of the seconds
    prec = 0;
end

sign_data = sign(data);
data = abs(data);

iDays = floor(data * 1 / 60 / 60 / 24);
iHours = floor(data * 1 / 60 / 60 - iDays * 24);
iMins = floor(data * 1 / 60 - iDays * 24 * 60 - iHours * 60 );
iSeconds = floor(data - iDays * 24 * 60  * 60 - iHours * 60  * 60 - iMins * 60);
iMilliSeconds = (data - iDays * 24 * 60  * 60 - iHours * 60  * 60 - iMins * 60 - iSeconds) * 1000;

ldata = length(data);
strRetVal = cell(ldata,1);

for ii = 1:ldata
    
    if( sign_data(ii) < 0 )
        strAux = '-'; 
    else
        strAux = [];
    end
    
    if( iDays(ii) > 0 )
        strAux = [ strAux num2str(iDays(ii)) 'd ' ]; 
    end
    
    if( iHours(ii) > 0 )
        strAux = [  strAux num2str(iHours(ii)) 'h ' ]; 
    end

    if( iMins(ii) > 0 )
        strAux = [  strAux num2str(iMins(ii)) ''' ' ]; 
    end

    if (iSeconds(ii) > 0  || isempty(strAux) )
        if( prec > 0 && iMilliSeconds(ii) > 0 )
            strAux = [  strAux sprintf( [ '%d.' ] , iSeconds(ii) ) ]; 
            strAux2 = sprintf( '%03d', round(iMilliSeconds(ii)) ); 
            strAux = [  strAux strAux2(1:min(3,prec)) '"' ]; 
        else
            strAux = [  strAux sprintf( [ '%d"' ] , iSeconds(ii)) ]; 
        end
    end
    
    strRetVal{ii} = strAux;
    
end

strRetVal = char(strRetVal);
