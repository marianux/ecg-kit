function strRetVal = Seconds2HMS(data, prec)

if( nargin < 2 )
    % decimals of the seconds
    prec = 0;
end

sign_data = sign(data);
data = abs(data);

iHours = floor(data * 1 / 60 / 60);
iMins = floor(data * 1 / 60 - iHours * 60 );
iSeconds = data - iHours * 60  * 60 - iMins * 60;
iMilliSeconds = iSeconds - floor(iSeconds);

if(iMilliSeconds == 0) 
    prec = 0;
end

ldata = length(data);
strRetVal = cell(ldata,1);

for ii = 1:ldata
    
    if( sign_data(ii) < 0 )
        strAux = '-'; 
    else
        strAux = [];
    end
    
    if( iHours(ii) > 0 )
        strAux = [   num2str(iHours(ii)) ' h ' ]; 
    end

    if( iMins(ii) > 0 )
        strAux = [  strAux num2str(iMins(ii)) ''' ' ]; 
    end

    if (iSeconds(ii) > 0  || isempty(strAux) )
        strAux = [  strAux sprintf( [ '%2.' num2str(prec) 'f"' ] , iSeconds(ii)) ]; 
    end
    
    strRetVal{ii} = strAux;
    
end

strRetVal = char(strRetVal);
