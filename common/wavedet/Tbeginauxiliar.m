% script which part of multilead T wave onset delineation
% Rute Almeida
% Last update: 07FEB2012 
%
% Designed for MATLAB Version R12; tested with MATLAB Version R13
Aon=[];
Von=[];
if ~isnan(picon)
    if picon<(intreg(i,1)+1) %#ok<IJCL>
        messages.warnings=[messages.warnings {'Ton out of the searching interval, in the fist step of multilead!'}];
        position.Ton(intervalo(i))=NaN;   %#ok<IJCL>
        position.Ttipoon(intervalo(i))=NaN;   %#ok<IJCL>
    else       
        if (pos.Ton-samp(1)+1)==picon
            begaux=(position0.Ton(intervalo(i))-samp(1)+1-round(messages.setup.wavedet.T_CSE_tol*messages.setup.wavedet.freq));   %#ok<IJCL>
        else
            begaux=max([intreg(i,1) (position0.Ton(intervalo(i))-samp(1)+1-messages.setup.wavedet.T_CSE_tol*messages.setup.wavedet.freq) ]);   %#ok<IJCL>
        end
        if ~isempty(w3)
            pontos2on=[w1(begaux:picon,scale) w2(begaux:picon,scale) w3(begaux:picon,scale)];
        else
            pontos2on=[w1(begaux:picon,scale) w2(begaux:picon,scale)];
        end
        weight2=ones(size(pontos2on,1),1);
        [aon,von] = optimline(pontos2on, weight2,OPT);
        
        Aon=[Aon;aon];
        Von=[Von;von];
        if ~isempty(w3)
            if i<length(timenew),   %#ok<IJCL>
                newleadbeatTon=(w1(timenew(i):timenew(i+1),scale).*von(1)+w2(timenew(i):timenew(i+1),scale)*von(2)+w3(timenew(i):timenew(i+1),scale)*von(3))./norm(von);   %#ok<IJCL>
               %signewTon=(sig(timenew(i):timenew(i+1),1).*von(1)+sig(timenew(i):timenew(i+1),2)*von(2)+sig(timenew(i):timenew(i+1),3)*von(3))./norm(von);
            else                        % if last beat of the segment
                newleadbeatTon=(w1(timenew(i):end,scale).*von(1)+w2(timenew(i):end,scale)*von(2)+w3(timenew(i):end,scale)*von(3))./norm(von);   %#ok<IJCL>
                %signewTon=(sig(timenew(i):end,1).*von(1)+sig(timenew(i):end,2)*von(2)+sig(timenew(i):end,3)*von(3))./norm(von);
            end
        else
            if i<length(timenew),  %#ok<IJCL>
                newleadbeatTon=(w1(timenew(i):timenew(i+1),scale).*von(1)+w2(timenew(i):timenew(i+1),scale)*von(2))./norm(von);   %#ok<IJCL>
                %signewTon=(sig(timenew(i):timenew(i+1),1).*von(1)+sig(timenew(i):timenew(i+1),2)*von(2))./norm(von);
            else                        % if last beat of the segment
                newleadbeatTon=(w1(timenew(i):end,scale).*von(1)+w2(timenew(i):end,scale)*von(2))./norm(von);   %#ok<IJCL>
                %signewTon=(sig(timenew(i):end,1).*von(1)+sig(timenew(i):end,2)*von(2))./norm(von);
            end
        end

        [posTon,piconTon2,picoffTon2,janelasTon2]=twave3D(newleadbeatTon,timenew,i,messages.setup.wavedet.freq,S,QRSoff,samp,rrmed,messages);   %#ok<IJCL>
        %janelasnew=janelas2+timenew(i)-1;
        piconall=[picon piconTon2];
        if isnan(piconTon2) | isempty(piconTon2) %#ok<OR2>
            amppiconall=[abs(newleadbeat(picon-timenew(i)+1)) NaN];   %#ok<IJCL>
        elseif isempty(piconTon2)
            amppiconall=[abs(newleadbeat(picon-timenew(i)+1)) NaN];   %#ok<IJCL>
        else
            amppiconall=[abs(newleadbeat(picon-timenew(i)+1)) abs(newleadbeatTon(piconTon2-timenew(i)+1))];   %#ok<IJCL>
        end
        newTon=[position0.Ton(intervalo(i)) posTon.Ton];   %#ok<IJCL>
        if (amppiconall(end)<amppiconall(end-1))
            position.contadorTon(intervalo(i))=-1;   %#ok<IJCL>
            position.Ton(intervalo(i))=position0.Ton(intervalo(i));   %#ok<IJCL>
            position.Ttipoon(intervalo(i))=position0.Ttipo(intervalo(i));   %#ok<IJCL>
        else
            position.contadorTon(intervalo(i))=0;   %#ok<IJCL>
            %%%%%%%%%%%%%%%%%%%%%%recursion
            while (~isnan(piconall(end)) &&  abs(newTon(end)-newTon(end-1))>Tconvergence_crit) && (amppiconall(end)>amppiconall(end-1))
                position.contadorTon(intervalo(i))=position.contadorTon(intervalo(i))+1;   %#ok<IJCL>
                if (pos.Ton-samp(1)+1)==picon
                    begaux=(pos.Ton-samp(1)+1-round(messages.setup.wavedet.T_CSE_tol*messages.setup.wavedet.freq));
                else
                    begaux=min([intreg(i,1) (pos.Ton-samp(1)+1-round(messages.setup.wavedet.T_CSE_tol*messages.setup.wavedet.freq))]);   %#ok<IJCL>
                end

                if ~isempty(w3)
                    pontosnewon=[w1(begaux:piconall(end),scale) w2(begaux:piconall(end),scale) w3(begaux:piconall(end),scale)];
                else
                    pontosnewon=[w1(begaux:piconall(end),scale) w2(begaux:piconall(end),scale)];
                end
                if size(pontosnewon,1)>1
                    weightnew=ones(size(pontosnewon,1),1);
                    [aon,von] = optimline(pontosnewon, weightnew,OPT);
                    Aon=[Aon;aon]; %#ok<AGROW>
                    Von=[Von;von]; %#ok<AGROW>

                    if ~isempty(w3)
                        if i<length(timenew),  %#ok<IJCL>
                            newleadbeatnewTon=(w1(timenew(i):timenew(i+1),scale).*von(1)+w2(timenew(i):timenew(i+1),scale)*von(2)+w3(timenew(i):timenew(i+1),scale)*von(3))./norm(von);  %#ok<IJCL>
                           % signewnewTon=(sig(timenew(i):timenew(i+1),1).*von(1)+sig(timenew(i):timenew(i+1),2)*von(2)+sig(timenew(i):timenew(i+1),3)*von(3))./norm(von);
                        else                        % if last beat of the segment
                            newleadbeatnewTon=(w1(timenew(i):end,scale).*von(1)+w2(timenew(i):end,scale)*von(2)+w3(timenew(i):end,scale)*von(3))./norm(von);   %#ok<IJCL>
                          %  signewnewTon=(sig(timenew(i):end,1).*von(1)+sig(timenew(i):end,2)*von(2)+sig(timenew(i):end,3)*von(3))./norm(von);
                        end
                    else
                        if i<length(timenew),   %#ok<IJCL>
                            newleadbeatnewTon=(w1(timenew(i):timenew(i+1),scale).*von(1)+w2(timenew(i):timenew(i+1),scale)*von(2))./norm(von);   %#ok<IJCL>
                           % signewnewTon=(sig(timenew(i):timenew(i+1),1).*von(1)+sig(timenew(i):timenew(i+1),2)*von(2))./norm(von);
                        else                        % if last beat of the segment
                            newleadbeatnewTon=(w1(timenew(i):end,scale).*von(1)+w2(timenew(i):end,scale)*von(2))./norm(von);   %#ok<IJCL>
                           % signewnewTon=(sig(timenew(i):end,1).*von(1)+sig(timenew(i):end,2)*von(2))./norm(von);
                        end
                    end
                    [posTon,piconTonnew,picoffTonnew,janelasTonnew,messages]=twave3D(newleadbeatnewTon,timenew,i,messages.setup.wavedet.freq,S,QRSoff,samp,rrmed,messages);   %#ok<IJCL>
                else
                    piconTonnew=NaN;
                end
                piconall=[piconall piconTonnew]; %#ok<AGROW>
                if isnan(piconTonnew)
                    amppiconall=[amppiconall NaN]; %#ok<AGROW>
                else
                    amppiconall=[amppiconall abs(newleadbeatnewTon(piconTonnew-timenew(i)+1))]; %#ok<AGROW,IJCL>
                end
                newTon=[newTon posTon.Ton]; %#ok<AGROW>

                if (amppiconall(end)<amppiconall(end-1))
                    posTon.Ton=NaN; % because the lead got worse
                end
                if isnan(posTon.Ton); %01Jun05
                    posTon.Ton=newTon(end-1);
                    position.contadorTon(intervalo(i))=position.contadorTon(intervalo(i))-1;   %#ok<IJCL>
                end
                if  prod(piconall(1:end-1)-piconall(end))==0  % ciclo infinito no teste 14
                    piconall(end)=NaN;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            position.Ton(intervalo(i))=posTon.Ton;   %#ok<IJCL>
            position.Ttipoon(intervalo(i))=posTon.Ttipo;   %#ok<IJCL>
                   
       end
        
    end

else
    position.Ton(intervalo(i))=NaN;   %#ok<IJCL>
    position.Ttipoon(intervalo(i))=NaN;   %#ok<IJCL>
end