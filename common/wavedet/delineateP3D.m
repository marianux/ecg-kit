function [position0,position,messages]= delineateP3D(heasig,w1,w2,w3,timenew,position0,position,intervalo,samp,ultimo_anot,messages)
%
% [position0,position,A,V,A2,V2]=% delineateP3D(heasig,w1,w2,w3,intreg,scale,timenew,position0,position,intervalo,samp,figuresoff,sig)
%
% Input Parameters:
%   heasig: struct vector with header information
%   w1: matrix with WT scales 1 to 5 for lead 1
%   w2: matrix with WT scales 1 to 5 for lead 2
%   w3: matrix with WT scales 1 to 5 for lead 1
%  intreg: interval to be considered in fitting the optimal direction
%  scale to use in delineation
%   timenew:  QRS times in inedexes refering the interval included in the current excerpt (borders excluded)
%   position0: struct vector with the detected points for 3D multilead  step1
%   position: struct vector with the detected points for 3D multilead  step2
%   position1: struct vector with the detected points for lead1 - to be replaced for multilead marks
%   samp: samples included in the current excerpt (borders excluded)
%   intervalo: numeration of the beats processed in this segment
%
%Output Parameters:
%   actualized parameters: position, position0
% A - points of fitted lines in step 1
% A2 - points of fitted lines in step 2
%  V- vectors director to the fitted line
%  V2- vectors director to the fitted line
%
% Rute Almeida: trabalho de doutoramento
% Last update: 07FEB2012 Rute Almeida
%
% MATLAB Version 6.5.0.180913a (R13)

if nargin<11
    messages = new_empty_mesg;
elseif nargin<10
    messages.errors=[messages.errors {'Fatal error in wavedet_3D\delineate3D: not enough inputs.'}];
    warning(char(messages.errors(end)))
    messages.errors_desc=[messages.errors_desc 'Mandatory inputs not defined.'];
    messages.status=0;
    return
end
if ~isfield(messages,'setup'), messages.setup.wavedet=[]; end
if ~isfield(messages,'errors'), messages.errors=[]; end
if ~isfield(messages,'errors_desc'), messages.errors_desc=[]; end
if ~isfield(messages,'warnings'), messages.warnings=[]; end
if isfield(messages,'status')
    if messages.status~=1
        messages.warnings=[messages.warnings {['Initial status=' num2str(messages.status) '. status changed to 1']}];
    end
end
messages.status=1;

global recursioncount OPT

%Threshols for onset and offset detection
%%%%%% Constants and Thresholds  !!!!!!!!!!!!!!!!!!!!!!!!!
if ~isfield(messages.setup.wavedet,'inivent_tol_P')
    messages.setup.wavedet.inivent_tol_P  =  0.34;
end
if ~isfield(messages.setup.wavedet,'finvent_tol_P')
    messages.setup.wavedet.finvent_tol_P  =  0.1;
end
if ~isfield(messages.setup.wavedet,'P_CSE_tol')
    messages.setup.wavedet.P_CSE_tol=10.2*2/1000; % based on CSE tolerance
end
if ~isfield(messages.setup.wavedet,'Pconvergence_crit')
    messages.setup.wavedet.Pconvergence_crit=0; % 1 maks are acepted as the same if they difer by Tconvergence_crit
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A2=[];
if ~isempty(w3)
    V=NaN*ones(size(timenew,2),3);
    A=NaN*ones(size(timenew,2),3);
else
    V=NaN*ones(size(timenew,2),2);
    A=NaN*ones(size(timenew,2),2);
end
V2=[];
intervalo=intervalo(1):intervalo(end);
for i=1:size(timenew,2)
    
    if ~isnan(position.QRSon(i+intervalo(1)-1))
        qrson = position.QRSon(i+intervalo(1)-1)-samp(1)+1;
    else
        qrson = timenew(i) - messages.setup.wavedet.finvent_tol_P*messages.setup.wavedet.freq;
    end
    
    ultimo_anot = ultimo_anot - samp(1)+1;
    % last anotation of the previous segment can overlap with the present one
    if (i==1 && ultimo_anot>=1),   % Protection in case it is the first beat
        begwin = max([1 qrson-round(messages.setup.wavedet.inivent_tol_P*messages.setup.wavedet.freq)+1 ultimo_anot+1]);
    else
        begwin = max(1, qrson-round(messages.setup.wavedet.inivent_tol_P*messages.setup.wavedet.freq)+1);
    end
    endwin = qrson- round(messages.setup.wavedet.finvent_tol_P*messages.setup.wavedet.freq)+1;
    
    if begwin>=endwin % empty search window: ultimo_anot should be too close to qrson
        position.Ptipoon(intervalo(i))=NaN;
        position.Ppicon(intervalo(i))=NaN;
        position.PscalePon(intervalo(i))=NaN;
        position.Pon(intervalo(i))=NaN;
        recursioncount.Pon(intervalo(i))=NaN;
        position.contadorPon(intervalo(i))=NaN;
        position.PscalePoff(intervalo(i))=NaN;
        position.Ppicoff(intervalo(i))=NaN;
        position.Poff(intervalo(i))=NaN;
        position.Ptipooff(intervalo(i))=NaN;
        recursioncount.Poff(intervalo(i))=NaN;
        position.contadorPoff(intervalo(i))=NaN; %Rute
        position.P(intervalo(i))=NaN;
    else
        
        scale=4;
        if ~isempty(w3)
            pontos=[w1(begwin:min(endwin,length(w1)),scale) w2(begwin:min(endwin,length(w2)),scale) w3(begwin:min(endwin,length(w3)),scale)];
        else
            pontos=[w1(begwin:min(endwin,length(w1)),scale) w2(begwin:min(endwin,length(w2)),scale)];
        end
        weight=ones(size(pontos,1),1);
        [a,v] = optimline(pontos, weight,OPT);
        A(i,:)=a;
        V(i,:)=v;
        %WT projected: from the previous QRS complex to the following QRS complex
        newleadbeat=[];
        if ~isempty(w3)
            if i>1,
                interval=timenew(i-1):timenew(i);
            else                        % if first beat of the segment
                interval=1:timenew(i);
            end
            
            newleadbeat(interval,3)=(w1(interval,3).*v(1)+w2(interval,3)*v(2)+w3(interval,3)*v(3))./norm(v); %#ok<AGROW>
            newleadbeat(interval,4)=(w1(interval,4).*v(1)+w2(interval,4)*v(2)+w3(interval,4)*v(3))./norm(v); %#ok<AGROW>
            newleadbeat(interval,5)=(w1(interval,5).*v(1)+w2(interval,5)*v(2)+w3(interval,5)*v(3))./norm(v); %#ok<AGROW>
            %signew(interval)=(sig(interval,1).*v(1)+sig(interval,2)*v(2)+sig(interval,3)*v(3))./norm(v);
        else
            if i>1,
                interval=timenew(i-1):timenew(i);
            else                        % if first beat of the segment
                interval=1:timenew(i);
            end
            % newleadbeat=NaN*ones(length(interval),5);
            newleadbeat(:,3)=(w1(interval,3).*v(1)+w2(timenew(i):timenew(i+1),3)*v(2))./norm(v);
            newleadbeat(:,4)=(w1(interval,4).*v(1)+w2(timenew(i):timenew(i+1),4)*v(2))./norm(v);
            newleadbeat(:,5)=(w1(interval,5).*v(1)+w2(timenew(i):timenew(i+1),5)*v(2))./norm(v);
            %signew=(sig(timenew(i):timenew(i+1),1).*v(1)+sig(timenew(i):timenew(i+1),2)*v(2))./norm(v);
        end
        if i>1,
            pos.qrs=position.qrs([i-1 i]);
            pos.QRSon=position.QRSon([i-1 i]);
            pos.QRSoff=position.QRSoff([i-1 i]);
            pos.Toff=position.Toff([i-1 i]);
            [pos,picon,picoff]= pwavef(heasig,samp,timenew([i i]),position,newleadbeat,[intervalo(i) intervalo(i)],ultimo_anot,messages);
        else
            pos.qrs=position.qrs(i);
            pos.QRSon=position.QRSon(i);
            pos.QRSoff=position.QRSoff(i);
            pos.Toff=position.Toff(i);
            [pos,picon,picoff]= pwavef(heasig,samp,timenew(i),position,newleadbeat,[intervalo(i) intervalo(i)] ,ultimo_anot,messages);
        end
        
        position0.Pscale(intervalo(i))=pos.Pscale(intervalo(i));
        position0.Pon(intervalo(i))=pos.Pon(intervalo(i));
        position0.Poff(intervalo(i))=pos.Poff(intervalo(i));
        position0.Ptipo(intervalo(i))=pos.Ptipo(intervalo(i));
        position0.direccionPonopt(:,intervalo(i))=NaN*ones(size(intervalo(i)));
        %position0.P(intervalo(i))=pos.P(intervalo(i));
        %position0.Pprima(intervalo(i))=pos.Pprima(intervalo(i));
        
        
        %%%%% P wave onset
        if ~isnan(picon)
            
            begaux=max([(position0.Pon(intervalo(i))-samp(1)+1-round(messages.setup.wavedet.P_CSE_tol*messages.setup.wavedet.freq)) 1]);
            if ~isempty(w3)
                pontos2on=[w1(begaux:picon,scale) w2(begaux:picon,scale) w3(begaux:picon,scale)];
            else
                pontos2on=[w1(begaux:picon,scale) w2(begaux:picon,scale)];
            end
            if size(pontos2on,1)<=1 % unable to do another step: empty search window
                position.contadorPon(intervalo(i))=-1;
                position.PscalePon(intervalo(i))=position0.Pscale(intervalo(i));
                position.Pon(intervalo(i))=position0.Pon(intervalo(i));
                position.Ppicon=picon+samp(1)-1;
                %position.P(intervalo(i))= position0.P(intervalo(i));
                %position.Pprima(intervalo(i))=position0.Pprima(intervalo(i));
                position.direccionPonopt(intervalo(i),:)=V(end,:)';
                position.Ptipoon(intervalo(i))=position0.Ptipo(intervalo(i));
            else
                weight2=ones(size(pontos2on,1),1);
                [a,v] = optimline(pontos2on, weight2,OPT);
                A2=[A2;a]; %#ok<AGROW>
                V2=[V2;v]; %#ok<AGROW>
                newleadbeat2on = [];
                if ~isempty(w3)
                    newleadbeat2on(interval,3)=(w1(interval,3).*v(1)+w2(interval,3)*v(2)+w3(interval,3)*v(3))./norm(v);%#ok<AGROW>
                    newleadbeat2on(interval,4)=(w1(interval,4).*v(1)+w2(interval,4)*v(2)+w3(interval,4)*v(3))./norm(v);%#ok<AGROW>
                    newleadbeat2on(interval,5)=(w1(interval,5).*v(1)+w2(interval,5)*v(2)+w3(interval,5)*v(3))./norm(v);%#ok<AGROW>
                    % signew2on(interval)=(sig(interval,1).*v(1)+sig(interval,2)*v(2)+sig(interval,3)*v(3))./norm(v);
                else
                    newleadbeat2on(interval,3)=(w1(interval,3).*v(1)+w2(timenew(i):timenew(i+1),3)*v(2))./norm(v);%#ok<AGROW>
                    newleadbeat2on(interval,4)=(w1(interval,4).*v(1)+w2(timenew(i):timenew(i+1),4)*v(2))./norm(v);%#ok<AGROW>
                    newleadbeat2on(interval,5)=(w1(interval,5).*v(1)+w2(timenew(i):timenew(i+1),5)*v(2))./norm(v);%#ok<AGROW>
                    %signew2on(interval)=(sig(timenew(i):timenew(i+1),1).*v(1)+sig(timenew(i):timenew(i+1),2)*v(2))./norm(v);
                end
                
                if i>1,
                    [pos2,picon2,picoff2]= pwavef(heasig,samp,timenew([i  i]),pos,newleadbeat2on,[intervalo(i) intervalo(i)],ultimo_anot,messages); %#ok<NASGU>
                else
                    [pos2,picon2,picoff2]= pwavef(heasig,samp,timenew(i),pos,newleadbeat2on,[intervalo(i) intervalo(i)] ,ultimo_anot,messages); %#ok<NASGU>
                end
                piconall=[picon picon2];
                if isnan(picon2)
                    amppiconall=[abs(newleadbeat(picon,scale)) NaN];
                else
                    amppiconall=[abs(newleadbeat(picon,scale)) abs(newleadbeat2on(picon2,scale))];
                end
                if (amppiconall(end)<amppiconall(end-1)) || isnan(amppiconall(end))% second step worse than first step
                    position.contadorPon(intervalo(i))=-1;
                    position.Ppicon(intervalo(i))=piconall(end-1)+samp(1)-1;
                    position.PscalePon(intervalo(i))=position0.Pscale(intervalo(i));
                    position.Pon(intervalo(i))=position0.Pon(intervalo(i));
                    position.direccionPonopt(:,intervalo(i))=V(end,:);
                    position.Ptipoon(intervalo(i))=position0.Ptipo(intervalo(i));
                elseif pos2.Pon(intervalo(i))==position0.Pon(intervalo(i)) % first and second steps give same mark
                    position.contadorPon(intervalo(i))=0;
                    position.Ppicon(intervalo(i))=piconall(end)+samp(1)-1;
                    position.PscalePon(intervalo(i))=pos2.Pscale(intervalo(i));
                    position.Pon(intervalo(i))=pos2.Pon(intervalo(i));
                    position.direccionPonopt(:,intervalo(i))=V2(end,:);
                    position.Ptipoon(intervalo(i))=pos2.Ptipo(intervalo(i));
                else % more steps
                    newPon=[position0.Pon(intervalo(i)) pos2.Pon(intervalo(i))];
                    position.contadorPon(intervalo(i))=0;
                    posnewon.Pscale(intervalo(i))=NaN;
                    posnewon.Pon(intervalo(i))=NaN;
                    posnewon.Ptipo(intervalo(i))=NaN;
                    position.Ppicon(intervalo(i))=piconall(end)+samp(1)-1;
                    while (~isnan(piconall(end)) &&  abs(newPon(end)-newPon(end-1))>messages.setup.wavedet.Pconvergence_crit) && (amppiconall(end)>amppiconall(end-1))
                        position.contadorPon(intervalo(i))=position.contadorPon(intervalo(i))+1;
                        begaux=max([(position0.Pon(intervalo(i))-samp(1)+1-round(messages.setup.wavedet.P_CSE_tol*messages.setup.wavedet.freq)) 1]);
                        if ~isempty(w3)
                            pontosnewon=[w1( begaux:piconall(end),scale) w2( begaux:piconall(end),scale) w3( begaux:piconall(end),scale)];
                        else
                            pontosnewon=[w1(begaux:piconall(end),scale) w2( begaux:piconall(end),scale)];
                        end
                        if size(pontosnewon,1)<=1 % unable to do another step: empty search window
                            posnewon.Pon(intervalo(i))=newPon(end-1); % previous step
                            posnewon.Pscale(intervalo(i))=pos2.Pscale(intervalo(i));
                            posnewon.Ptipo(intervalo(i))=pos2.Ptipo(intervalo(i));
                            position.contadorPon(intervalo(i))=position.contadorPon(intervalo(i))-1;
                            piconall=[piconall NaN]; %#ok<AGROW>
                            amppiconall=[amppiconall NaN]; %#ok<AGROW>
                        else
                            weight2=ones(size(pontosnewon,1),1);
                            [a,v] = optimline(pontosnewon, weight2,OPT);
                            A2=[A2;a]; %#ok<AGROW>
                            V2=[V2;v]; %#ok<AGROW>
                            newleadbeatnewon = [];
                            if ~isempty(w3)
                                % newleadbeatnewon=NaN*ones(interval(end),5);
                                newleadbeatnewon(interval,3)=(w1(interval,3).*v(1)+w2(interval,3)*v(2)+w3(interval,3)*v(3))./norm(v);%#ok<AGROW>
                                newleadbeatnewon(interval,4)=(w1(interval,4).*v(1)+w2(interval,4)*v(2)+w3(interval,4)*v(3))./norm(v);%#ok<AGROW>
                                newleadbeatnewon(interval,5)=(w1(interval,5).*v(1)+w2(interval,5)*v(2)+w3(interval,5)*v(3))./norm(v);%#ok<AGROW>
                                %signewon(interval)=(sig(interval,1).*v(1)+sig(interval,2)*v(2)+sig(interval,3)*v(3))./norm(v);
                            else
                                % newleadbeatnewon=NaN*ones(interval(end),5);
                                newleadbeatnewon(interval,3)=(w1(interval,3).*v(1)+w2(timenew(i):timenew(i+1),3)*v(2))./norm(v);%#ok<AGROW>
                                newleadbeatnewon(interval,4)=(w1(interval,4).*v(1)+w2(timenew(i):timenew(i+1),4)*v(2))./norm(v);%#ok<AGROW>
                                newleadbeatnewon(interval,5)=(w1(interval,5).*v(1)+w2(timenew(i):timenew(i+1),5)*v(2))./norm(v);%#ok<AGROW>
                                %signewon(interval)=(sig(timenew(i):timenew(i+1),1).*v(1)+sig(timenew(i):timenew(i+1),2)*v(2))./norm(v);
                            end
                            if i>1,
                                [posnewon,piconnew,picoffnew]= pwavef(heasig,samp,timenew([i  i]),pos,newleadbeatnewon,[intervalo(i) intervalo(i)],ultimo_anot,messages); %#ok<NASGU>
                            else
                                [posnewon,piconnew,picoffnew]= pwavef(heasig,samp,timenew(i),pos,newleadbeatnewon,[intervalo(i) intervalo(i)] ,ultimo_anot,messages); %#ok<NASGU>
                            end
                            newPon=[newPon posnewon.Pon(intervalo(i))]; %#ok<AGROW>
                            if isnan(piconnew)
                                amppiconall=[amppiconall NaN]; %#ok<AGROW>
                            else
                                amppiconall=[amppiconall abs(newleadbeatnewon(piconnew,scale))]; %#ok<AGROW>
                            end
                            piconall=[piconall piconnew]; %#ok<AGROW>
                            
                            if isnan(posnewon.Pon(intervalo(i))) || (amppiconall(end)<amppiconall(end-1)) % lead got worse
                                posnewon.Pon(intervalo(i))=newPon(end-1); % previous step
                                position.Ppicon(intervalo(i))=piconall(end-1)+samp(1)-1;
                                posnewon.Pscale(intervalo(i))=pos2.Pscale(intervalo(i));
                                posnewon.Ptipo(intervalo(i))=pos2.Ptipo(intervalo(i));
                                position.contadorPon(intervalo(i))=position.contadorPon(intervalo(i))-1;
                            else
                                position.direccionPonopt(:,intervalo(i))=V2(end,:);
                                position.Ppicon(intervalo(i))=piconall(end)+samp(1)-1;
                            end
                            if  prod(piconall(1:end-1)-piconall(end))==0  % ciclo infinito no teste 14
                                piconall(end)=NaN;
                            end
                            
                        end % if size(pontosnewon,1)>1
                    end %while
                    position.PscalePon(intervalo(i))=posnewon.Pscale(intervalo(i));
                    position.Pon(intervalo(i))=posnewon.Pon(intervalo(i));
                    position.Ptipoon(intervalo(i))=posnewon.Ptipo(intervalo(i));
                end % more steps
            end   %size(pontosnewon,1)>1
        else % isnan(picon)
            position.Ptipoon(intervalo(i))=NaN;
            position.Ppicon(intervalo(i))=NaN;
            position.PscalePon(intervalo(i))=NaN;
            position.Pon(intervalo(i))=NaN;
            recursioncount.Pon(intervalo(i))=NaN;
            position.contadorPon(intervalo(i))=NaN; %Rute
        end
        %%%%%% P wave onset
        %%%%% P wave end
        if ~isnan(picoff)
            
            endaux=min([(position0.Poff(intervalo(i))-samp(1)+1+round(messages.setup.wavedet.P_CSE_tol*messages.setup.wavedet.freq)) length(w1)]);
            if ~isempty(w3)
                pontos2off=[w1(picoff:endaux,scale) w2(picoff:endaux,scale) w3(picoff:endaux,scale)];
            else
                pontos2off=[w1(picoff:endaux,scale) w2(picoff:endaux,scale)];
            end
            if size(pontos2off,1)<=1 % unable to do another step: empty search window
                position.contadorPoff(intervalo(i))=-1;
                position.PscalePoff(intervalo(i))=position0.Pscale(intervalo(i));
                position.Poff(intervalo(i))=position0.Poff(intervalo(i));
                position.Ppicoff=picoff+samp(1)-1;
                %position.P(intervalo(i))= position0.P(intervalo(i));
                %position.Pprima(intervalo(i))=position0.Pprima(intervalo(i));
                position.direccionPoffopt(:,intervalo(i))=V(end,:);
                position.Ptipooff(intervalo(i))=position0.Ptipo(intervalo(i));
            else
                weight2=ones(size(pontos2off,1),1);
                [a,v] = optimline(pontos2off, weight2,OPT);
                A2=[A2;a]; %#ok<AGROW>
                V2=[V2;v]; %#ok<AGROW>
                newleadbeat2off = [];
                if ~isempty(w3)
                    %   newleadbeat2off=NaN*ones(interval(end),5);
                    newleadbeat2off(interval,3)=(w1(interval,3).*v(1)+w2(interval,3)*v(2)+w3(interval,3)*v(3))./norm(v);%#ok<AGROW>
                    newleadbeat2off(interval,4)=(w1(interval,4).*v(1)+w2(interval,4)*v(2)+w3(interval,4)*v(3))./norm(v);%#ok<AGROW>
                    newleadbeat2off(interval,5)=(w1(interval,5).*v(1)+w2(interval,5)*v(2)+w3(interval,5)*v(3))./norm(v);%#ok<AGROW>
                    %signew2off(interval)=(sig(interval,1).*v(1)+sig(interval,2)*v(2)+sig(interval,3)*v(3))./norm(v);
                else
                    % newleadbeat2off=NaN*ones(interval(end),5);
                    newleadbeat2off(interval,3)=(w1(interval,3).*v(1)+w2(timenew(i):timenew(i+1),3)*v(2))./norm(v);%#ok<AGROW>
                    newleadbeat2off(interval,4)=(w1(interval,4).*v(1)+w2(timenew(i):timenew(i+1),4)*v(2))./norm(v);%#ok<AGROW>
                    newleadbeat2off(interval,5)=(w1(interval,5).*v(1)+w2(timenew(i):timenew(i+1),5)*v(2))./norm(v);%#ok<AGROW>
                    %signew2off(interval)=(sig(timenew(i):timenew(i+1),1).*v(1)+sig(timenew(i):timenew(i+1),2)*v(2))./norm(v);
                end
                
                if i>1,
                    [pos2off,picon2off,picoff2]= pwavef(heasig,samp,timenew([i  i]),pos,newleadbeat2off,[intervalo(i) intervalo(i)],ultimo_anot,messages); %#ok<ASGLU>
                else
                    [pos2off,picon2off,picoff2]= pwavef(heasig,samp,timenew(i),pos,newleadbeat2off,[intervalo(i) intervalo(i)] ,ultimo_anot,messages); %#ok<ASGLU>
                end
                
                picoffall=[picoff picoff2];
                if isnan(picoff2)
                    amppicoffall=[abs(newleadbeat(picoff,scale)) NaN];
                else
                    amppicoffall=[abs(newleadbeat(picoff,scale)) abs(newleadbeat2off(picoff2,scale))];
                end
                if (amppicoffall(end)<amppicoffall(end-1))  || isnan(amppicoffall(end)) % second step worse than first step
                    position.contadorPoff(intervalo(i))=-1;
                    position.Ppicoff(intervalo(i))=picoffall(end-1)+samp(1)-1;
                    position.PscalePoff(intervalo(i))=position0.Pscale(intervalo(i));
                    position.Poff(intervalo(i))=position0.Poff(intervalo(i));
                    %position.P(intervalo(i))= position0.P(intervalo(i));
                    %position.Pprima(intervalo(i))=position0.Pprima(intervalo(i));
                    position.direccionPoffopt(:,intervalo(i))=V(end,:);
                    position.Ptipooff(intervalo(i))=position0.Ptipo(intervalo(i));
                elseif pos2off.Poff(intervalo(i))==position0.Poff(intervalo(i)) % first and second steps give same mark
                    position.contadorPoff(intervalo(i))=0;
                    position.Ppicoff(intervalo(i))=picoffall(end)+samp(1)-1;
                    position.PscalePoff(intervalo(i))=pos2off.Pscale(intervalo(i));
                    position.Poff(intervalo(i))=pos2off.Poff(intervalo(i));
                    %position.P(intervalo(i))= pos2off.P(intervalo(i));
                    %position.Pprima(intervalo(i))=pos2off.Pprima(intervalo(i));
                    position.direccionPoffopt(:,intervalo(i))=V2(end,:);
                    position.Ptipooff(intervalo(i))=pos2off.Ptipo(intervalo(i));
                else % more steps
                    newPoff=[position0.Poff(intervalo(i)) pos2off.Poff(intervalo(i))];
                    position.contadorPoff(intervalo(i))=0;
                    position.Ppicoff(intervalo(i))=picoffall(end)+samp(1)-1;
                    posnewoff.Pscale(intervalo(i))=NaN;
                    posnewoff.Poff(intervalo(i))=NaN;
                    posnewoff.Ptipo(intervalo(i))=NaN;
                    while (~isnan(picoffall(end)) &&  abs(newPoff(end)-newPoff(end-1))>messages.setup.wavedet.Pconvergence_crit ) && (amppicoffall(end)>amppicoffall(end-1))
                        position.contadorPoff(intervalo(i))=position.contadorPoff(intervalo(i))+1;
                        endaux=min([(position0.Poff(intervalo(i))-samp(1)+1+round(messages.setup.wavedet.P_CSE_tol*messages.setup.wavedet.freq)) length(w1)]);
                        if ~isempty(w3)
                            pontosnewoff=[w1(picoffall(end):endaux,scale) w2(picoffall(end):endaux,scale) w3(picoffall(end):endaux,scale)];
                        else
                            pontosnewoff=[w1(picoffall(end):endaux,scale) w2(picoffall(end):endaux,scale)];
                        end
                        if size(pontosnewoff,1)<=1 % unable to do another step: empty search window
                            posnewoff.Poff(intervalo(i))=newPoff(end-1); % previous step
                            posnewoff.Pscale(intervalo(i))=pos2off.Pscale(intervalo(i));
                            posnewoff.Ptipo(intervalo(i))=pos2off.Ptipo(intervalo(i));
                            %posnew.P=pos2off.P;
                            %posnew.Pprima=pos2off.Pprima;
                            position.contadorPoff(intervalo(i))=position.contadorPoff(intervalo(i))-1;
                            picoffall=[picoffall NaN]; %#ok<AGROW>
                            amppicoffall=[amppicoffall NaN]; %#ok<AGROW>
                        else
                            weight2=ones(size(pontosnewoff,1),1);
                            [a,v] = optimline(pontosnewoff, weight2,OPT);
                            A2=[A2;a]; %#ok<AGROW>
                            V2=[V2;v]; %#ok<AGROW>
                            newleadbeatnewoff = [];
                            if ~isempty(w3)
                                %newleadbeatnewoff=NaN*ones(interval(end),5);
                                newleadbeatnewoff(interval,3)=(w1(interval,3).*v(1)+w2(interval,3)*v(2)+w3(interval,3)*v(3))./norm(v);%#ok<AGROW>
                                newleadbeatnewoff(interval,4)=(w1(interval,4).*v(1)+w2(interval,4)*v(2)+w3(interval,4)*v(3))./norm(v);%#ok<AGROW>
                                newleadbeatnewoff(interval,5)=(w1(interval,5).*v(1)+w2(interval,5)*v(2)+w3(interval,5)*v(3))./norm(v);%#ok<AGROW>
                                %signewoff(interval)=(sig(interval,1).*v(1)+sig(interval,2)*v(2)+sig(interval,3)*v(3))./norm(v);
                            else
                                %newleadbeatnewoff=NaN*ones(interval(end),5);
                                newleadbeatnewoff(interval,3)=(w1(interval,3).*v(1)+w2(timenew(i):timenew(i+1),3)*v(2))./norm(v);%#ok<AGROW>
                                newleadbeatnewoff(interval,4)=(w1(interval,4).*v(1)+w2(timenew(i):timenew(i+1),4)*v(2))./norm(v);%#ok<AGROW>
                                newleadbeatnewoff(interval,5)=(w1(interval,5).*v(1)+w2(timenew(i):timenew(i+1),5)*v(2))./norm(v);%#ok<AGROW>
                                %signewoff(interval)=(sig(timenew(i):timenew(i+1),1).*v(1)+sig(timenew(i):timenew(i+1),2)*v(2))./norm(v);
                            end
                            if i>1,
                                [posnewoff,piconnew,picoffnew]= pwavef(heasig,samp,timenew([i  i]),pos,newleadbeatnewoff,[intervalo(i) intervalo(i)],ultimo_anot,messages); %#ok<ASGLU>
                            else
                                [posnewoff,piconnew,picoffnew]= pwavef(heasig,samp,timenew(i),pos,newleadbeatnewoff,[intervalo(i) intervalo(i)] ,ultimo_anot,messages); %#ok<ASGLU>
                            end
                            newPoff=[newPoff posnewoff.Poff(intervalo(i))]; %#ok<AGROW>
                            if isnan(picoffnew)
                                amppicoffall=[amppicoffall NaN]; %#ok<AGROW>
                            else
                                amppicoffall=[amppicoffall abs(newleadbeatnewoff(picoffnew,scale))]; %#ok<AGROW>
                            end
                            picoffall=[picoffall picoffnew]; %#ok<AGROW>
                            
                            
                            if isnan(posnewoff.Poff(intervalo(i))) || (amppicoffall(end)<amppicoffall(end-1)) % lead got worse
                                position.Ppicoff(intervalo(i))=picoffall(end-1)+samp(1)-1;
                                posnewoff.Poff(intervalo(i))=newPoff(end-1); % previous step
                                posnewoff.Pscale(intervalo(i))=pos2off.Pscale(intervalo(i));
                                posnewoff.Ptipo(intervalo(i))=pos2off.Ptipo(intervalo(i));
                                %posnewoff.P=pos2.P;
                                %posnewoff.Pprima=pos2.Pprima;
                                position.contadorPoff(intervalo(i))=position.contadorPoff(intervalo(i))-1;
                            else
                                position.direccionPoffopt(:,intervalo(i))=V2(end,:);
                                position.Ppicoff(intervalo(i))=picoffall(end)+samp(1)-1;
                            end
                            
                            if  prod(picoffall(1:end-1)-picoffall(end))==0
                                picoffall(end)=NaN;
                            end
                        end % if size(pontosnewoff,1)>1
                    end %while
                    position.PscalePoff(intervalo(i))=posnewoff.Pscale(intervalo(i));
                    position.Poff(intervalo(i))=posnewoff.Poff(intervalo(i));
                    position.Ptipooff(intervalo(i))=posnewoff.Ptipo(intervalo(i));
                    % position.P(intervalo(i))= posnewoff.P(intervalo(i));
                    % position.Pprima(intervalo(i))=posnewoff.Pprima(intervalo(i));
                end % more steps
            end   %size(pontosnewoff,1)>1
        else % isnan(picoff)
            position.PscalePoff(intervalo(i))=NaN;
            position.Ppicoff(intervalo(i))=NaN;
            position.Poff(intervalo(i))=NaN;
            position.Ptipooff(intervalo(i))=NaN;
            %position.P(intervalo(i))=NaN;
            %position.Pprima(intervalo(i))=NaN;
            recursioncount.Poff(intervalo(i))=NaN;
            position.contadorPoff(intervalo(i))=NaN; %Rute
        end
        %%%%%% P wave end
        
        %P peak
        if ~isnan(position.Ppicon(intervalo(i))) && ~isnan(position.Ppicoff(intervalo(i))) && position.Ppicon(intervalo(i))<position.Ppicoff(intervalo(i))
            pontosP=[w1((position.Ppicon(intervalo(i))-samp(1)+1):(position.Ppicoff(intervalo(i))-samp(1)+1),scale) w2((position.Ppicon(intervalo(i))-samp(1)+1):(position.Ppicoff(intervalo(i))-samp(1)+1),scale) w3((position.Ppicon(intervalo(i))-samp(1)+1):(position.Ppicoff(intervalo(i))-samp(1)+1),scale)];
            [m1,m2]=min((pontosP(:,1)).^2+(pontosP(:,2)).^2+(pontosP(:,3)).^2); %#ok<ASGLU> % main P peak
            position.P(intervalo(i))=m2+position.Ppicon(intervalo(i));
        else
            position.P(intervalo(i))=position0.P(intervalo(i));
        end
        %%
        %writennot protection: at least one sample between peak and border
        if (position.P(intervalo(i))-position.Pon(intervalo(i)))<1
            position.P(intervalo(i))=NaN;
        end
        if (position.Poff(intervalo(i))-position.P(intervalo(i)))<1
            position.P(intervalo(i))=NaN;
        end
        clear pos pos2 pontos pontos2on pontos2off posnew pos2off posnewoff
    end
end % for
