function [position0,position,messages]= multiqrsend(heasig,w1,w2,w3,auxqrspicoffset,QRSoff,time,timenew,position0,position,intervalo,samp,messages)
%[position0,position]= multiqrsend(heasig,w1,w2,w3,auxqrspiconset,QRSon,time,timenew,position0,position,intervalo,position1,position2,position3,samp)
%
% INPUT:
%   heasig: header information
%   w1,w2,w3: WT scales 1 to 3
%   auxqrspicoffset: latest relevant modulus maximum location associated to QRS complex found in any lead
%   QRSoff: latest QRS end location in any lead (initial QRS end location)
%   time: beats processed in this segment in each lead
%   timenew: beats processed in this segment (multilead)
%   position: struct vector with the detected points at the first multilead
%   step
%   position: struct vector with the detected points
%   intervalo: numeration of the beats processed in this segment
%   position1,position2,position3: struct vectors with the detected points using each specified lead
%   samp: samples included in the current excerpt (borders excluded)
%
% OUTPUT:
% actualized position0,position
%
% Rute Almeida  
% Last update: 07FEB2012
global recursioncount OPT
if ~isfield(messages.setup.wavedet,'Kron')
    messages.setup.wavedet.Kron = 20;%      3.5 5 8 10 15 20!
end
if ~isfield(messages.setup.wavedet,'Kroff')
    messages.setup.wavedet.Kroff = 14;
end
if ~isfield(messages.setup.wavedet,'Kq')% Q onset, antiguo Kq==1;                2 2.5 3.5 5 10 15!
    messages.setup.wavedet.Kq = 15;
end
if ~isfield(messages.setup.wavedet,'Ks')
    messages.setup.wavedet.Ks = 8;	      % S offset (Proyeuto de Mapi) .          2.5 3 5 8;
end
if ~isfield(messages.setup.wavedet,'maxQRSoff')
    messages.setup.wavedet.maxQRSoff =11.6*2; % based on CSE tolerance
end
% if ~isfield(messages.setup.wavedet,'tolsamepeak')
%         messages.setup.wavedet.tolsamepeak=5; %%%% number in samples to consider that a different peak was chosen      
% end
% tolsamepeak=messages.setup.wavedet.tolsamepeak;

% Kron=messages.setup.wavedet.Kron;
% Kroff=messages.setup.wavedet.Kroff;
% Kq=messages.setup.wavedet.Kq;
% Ks=messages.setup.wavedet.Ks;
maxQRSoff=messages.setup.wavedet.maxQRSoff;

auxqrsoff=nanmax(QRSoff-samp(1)+1)+round(maxQRSoff/1000*messages.setup.wavedet.freq);
for i=1:size(time,2)
    position.QRSoffcriteria((i+intervalo(1)-1))=2;% default % 1 not considered
    scale=2;
    if size(w1,2)==3
        pontos=[w1(auxqrspicoffset(i):min(auxqrsoff(i),end),scale)... 
            w2(auxqrspicoffset(i):min(auxqrsoff(i),end),scale)...
            w3(auxqrspicoffset(i):min(auxqrsoff(i),end),scale)]; %4JUL2011 % protect with end
    else
        pontos=[w1(auxqrspicoffset(i):min(auxqrsoff(i),end),scale)...
            w2(auxqrspicoffset(i):min(auxqrsoff(i),end),scale)];%4JUL2011 % protect with end
    end
    weight=ones(size(pontos,1),1);
    [a,v] = optimline(pontos, weight,OPT); %#ok<ASGLU>
    time0=round(nanmedian(time));% Rute 01.03.2010
    waux=nan*ones(size(w1)); %Rute6NOV08
    if i==1
        if size(timenew,2)>1 %4JUL08
            intervalonewlead=1:timenew(i+1);
        else
            intervalonewlead=1:length(w1);
        end
    elseif i<length(timenew)
        intervalonewlead=timenew(i-1):timenew(i+1);
    else
        intervalonewlead=timenew(i-1):length(w1);
    end
    if size(w1,2)==3
        newleadbeat1=(w1(intervalonewlead,1).*v(1)+w2(intervalonewlead,1)*v(2)+w3(intervalonewlead,1)*v(3))./norm(v);
        newleadbeat2=(w1(intervalonewlead,2).*v(1)+w2(intervalonewlead,2)*v(2)+w3(intervalonewlead,2)*v(3))./norm(v);
        %signew=(sig(intervalonewlead,1).*v(1)+sig(intervalonewlead,2)*v(2)+sig(intervalonewlead,3)*v(3))./norm(v);
    else
        newleadbeat1=(w1(intervalonewlead,1).*v(1)+w2(intervalonewlead,1)*v(2))./norm(v);
        newleadbeat2=(w1(intervalonewlead,2).*v(1)+w2(intervalonewlead,2)*v(2))./norm(v);
        %signew=(sig(intervalonewlead,1).*v(1)+sig(intervalonewlead,2)*v(2))./norm(v);
    end
    waux(intervalonewlead,1)=newleadbeat1;
    waux(intervalonewlead,2)=newleadbeat2;
    [position0,qrspiconset0,qrspicoffset0]= qrswavef(heasig,samp,time0(i),position0,waux,[i+intervalo(1)-1 i+intervalo(1)-1],messages); %#ok<ASGLU>
    if isempty(qrspicoffset0)||isnan(qrspicoffset0) % RUTE 6ABRI09
%         picoffall= NaN; %#ok<NASGU>
%         amppicoffall= NaN; %#ok<NASGU>
%         amppicoffaux= NaN; %#ok<NASGU>
        recursioncount.QRSoff((i+intervalo(1)-1))=-2;
        position.QRSoff(i+intervalo(1)-1) = NaN;
        position.S(i+intervalo(1)-1) = NaN;
        position_safe.R(i+intervalo(1)-1) = NaN; %13.AGO.07
        position.RR_inQRSoff(i+intervalo(1)-1) = NaN;
        position.Rprima(i+intervalo(1)-1)=NaN;
         
    else
        if size(w1,2)==3
            pontos2=[w1(qrspicoffset0:min(end,(position0.QRSoff(i+intervalo(1)-1)-samp(1)+1+round(maxQRSoff/1000*messages.setup.wavedet.freq))),scale)...
                w2(qrspicoffset0:min(end,(position0.QRSoff(i+intervalo(1)-1)-samp(1)+1+round(maxQRSoff/1000*messages.setup.wavedet.freq))),scale)...
                w3(qrspicoffset0:min(end,(position0.QRSoff(i+intervalo(1)-1)-samp(1)+1+round(maxQRSoff/1000*messages.setup.wavedet.freq))),scale)];
        else
            pontos2=[w1(qrspicoffset0:min(end,(position0.QRSoff(i+intervalo(1)-1)-samp(1)+1+round(maxQRSoff/1000*messages.setup.wavedet.freq))),scale)...
                w2(qrspicoffset0:min(end,(position0.QRSoff(i+intervalo(1)-1)-samp(1)+1+round(maxQRSoff/1000*messages.setup.wavedet.freq))),scale)];
        end
        weight=ones(length(pontos2),1);
        [a2,v2] = optimline(pontos2, weight,OPT);        %#ok<ASGLU>
        if i==1
            if size(timenew,2)>1 %4JUL08
                intervalonewlead=1:timenew(i+1);
            else
                intervalonewlead=1:length(w1);
            end
        elseif i<length(timenew)
            intervalonewlead=timenew(i-1):timenew(i+1);
        else
            intervalonewlead=timenew(i-1):length(w1);
        end
        if size(w1,2)==3
            newleadbeat12=(w1(intervalonewlead,1).*v2(1)+w2(intervalonewlead,1)*v2(2)+w3(intervalonewlead,1)*v2(3))./norm(v2);
            newleadbeat22=(w1(intervalonewlead,2).*v2(1)+w2(intervalonewlead,2)*v2(2)+w3(intervalonewlead,2)*v2(3))./norm(v2);
            %signew2=(sig(intervalonewlead,1).*v2(1)+sig(intervalonewlead,2)*v2(2)+sig(intervalonewlead,3)*v2(3))./norm(v2);
        else
            newleadbeat12=(w1(intervalonewlead,1).*v2(1)+w2(intervalonewlead,1)*v2(2))./norm(v2);
            newleadbeat22=(w1(intervalonewlead,2).*v2(1)+w2(intervalonewlead,2)*v2(2))./norm(v2);
            %signew2=(sig(intervalonewlead,1).*v2(1)+sig(intervalonewlead,2)*v2(2))./norm(v2);
        end
        waux=nan*ones(size(w1)); %Rute6NOV08
        waux(intervalonewlead,1)=newleadbeat12;
        waux(intervalonewlead,2)=newleadbeat22;

        [positionaux,qrspiconset2,qrspicoffset2]= qrswavef(heasig,samp,time0(i),position,waux,[i+intervalo(1)-1 i+intervalo(1)-1],messages); %#ok<ASGLU>
        position.QRSoff(i+intervalo(1)-1) = positionaux.QRSoff(i+intervalo(1)-1);
        position.S(i+intervalo(1)-1) = positionaux.S(i+intervalo(1)-1);
        position.R_inQRSoff(i+intervalo(1)-1) = positionaux.R(i+intervalo(1)-1) ;
        position_safe.R(i+intervalo(1)-1) = positionaux.R(i+intervalo(1)-1) ;%13.AGO.07
        position.Rprima(i+intervalo(1)-1) = positionaux.Rprima(i+intervalo(1)-1);

        if isempty(qrspicoffset2)||isnan(qrspicoffset2) %16MAY2011
            picoffall=[qrspicoffset0 NaN];
            amppicoffall=[(newleadbeat2(qrspicoffset0-intervalonewlead(1)+1)) NaN];
            amppicoffaux= NaN;
        else
            picoffall=[qrspicoffset0 qrspicoffset2];
            amppicoffall=[(newleadbeat2(qrspicoffset0-intervalonewlead(1)+1)) (newleadbeat22(qrspicoffset2-intervalonewlead(1)+1))];
            amppicoffaux= newleadbeat2(qrspicoffset2-intervalonewlead(1)+1);
        end
        newQRSoff=[NaN position0.QRSoff(i+intervalo(1)-1) position.QRSoff(i+intervalo(1)-1)];
        newleadbeat22aux=newleadbeat22;
        % stop criteria iii)  |  iv) ou novo
        %  if (positionaux.QRSoff(i+intervalo(1)-1)<=position.R(i+intervalo(1)-1)&...
        %position0.QRSoff(i+intervalo(1)-1)>position.R(i+intervalo(1)-1))| ...
        %((abs(amppicoffall(end))<abs(amppicoffall(end-1)))&((amppicoffaux(end)*amppicoffall(end))>0)  &...
        %~isnan(position0.S(i+intervalo(1)-1)))|( isnan(picoffall(end))) %| (abs(newQRSoff(end)-newQRSoff(end-1))<=1) % RUTE 9.07.2007
        %13.AGO.07
        % stop criteria iii)  |  iv) % opçao0 % parte da alternativa2
        if ((abs(amppicoffall(end))<abs(amppicoffall(end-1)))&&((amppicoffaux(end)*amppicoffall(end))>0)  &&...
                ~isnan(position0.S(i+intervalo(1)-1)))||( isnan(picoffall(end))) %| (abs(newQRSoff(end)-newQRSoff(end-1))<=1) % RUTE 9.07.2007
            % stop criteria iii)  | iv) | i) % alternativa1
            %if ((abs(amppicoffall(end))<abs(amppicoffall(end-1)))&((amppicoffaux(end)*amppicoffall(end))>0)  &...
            %~isnan(position0.S(i+intervalo(1)-1)))|( isnan(picoffall(end))) | (abs(newQRSoff(end)-newQRSoff(end-1))<=1) % RUTE 9.07.2007
            %alternative criteria stop criteria iii)without the S wave condition| iv) alternativa3
            %if ((abs(amppicoffall(end))<abs(amppicoffall(end-1)))&((amppicoffaux(end)*amppicoffall(end))>0) )|( isnan(picoffall(end))) %| (abs(newQRSoff(end)-newQRSoff(end-1))<=1) % RUTE 9.07.2007
            %if ( isnan(picoffall(end))) | (abs(newQRSoff(end)-newQRSoff(end-1))<=1) % RUTE 9.07.2007
            if isnan(picoffall(end))
                position.QRSoffcriteria(i+intervalo(1)-1)=4;
                % alternativa2 - exclude criteria i)
            elseif (abs(newQRSoff(end)-newQRSoff(end-1))<=1)
                position.QRSoffcriteria(i+intervalo(1)-1)=1;
            else
                position.QRSoffcriteria(i+intervalo(1)-1)=3;
            end
            recursioncount.QRSoff((i+intervalo(1)-1))=-1;
            position.QRSoff(i+intervalo(1)-1) = position0.QRSoff(i+intervalo(1)-1);
            position.S(i+intervalo(1)-1) = position0.S(i+intervalo(1)-1);
            position_safe.R(i+intervalo(1)-1) = position0.R(i+intervalo(1)-1) ; %13.AGO.07
            position.Rprima(i+intervalo(1)-1) = position0.Rprima(i+intervalo(1)-1);
        else
            recursioncount.QRSoff((i+intervalo(1)-1))=0;
            % stop criteria  iv)  | ii) | i) % opçao0 % parte da alternativa1
            %while ~isnan(picoffall(end)) &  nansum(newQRSoff(end)-newQRSoff(1:end-1)==0)<2 & (abs(newQRSoff(end)-newQRSoff(end-1))>1)%%%
            %alternative criteria stop criteria iv)  | ii)) % alternativa2e3
            while ~isnan(picoffall(end)) &&  nansum(newQRSoff(end)-newQRSoff(1:end-1)==0)<2 %& (abs(newQRSoff(end)-newQRSoff(end-1))>1)%%%
                recursioncount.QRSoff((i+intervalo(1)-1))=recursioncount.QRSoff((i+intervalo(1)-1))+1;
                if size(w1,2)==3
                    pontosnew=[w1(qrspicoffset2:min(end,(position.QRSoff(i+intervalo(1)-1)-samp(1)+1+round(maxQRSoff/1000*messages.setup.wavedet.freq))),scale)...
                        w2(qrspicoffset2:min(end,(position.QRSoff(i+intervalo(1)-1)-samp(1)+1+round(maxQRSoff/1000*messages.setup.wavedet.freq))),scale)...
                        w3(qrspicoffset2:min(end,(position.QRSoff(i+intervalo(1)-1)-samp(1)+1+round(maxQRSoff/1000*messages.setup.wavedet.freq))),scale)];
                else
                    pontosnew=[w1(qrspicoffset2:min(end,(position.QRSoff(i+intervalo(1)-1)-samp(1)+1+round(maxQRSoff/1000*messages.setup.wavedet.freq))),scale)...
                        w2(qrspicoffset2:min(end,(position.QRSoff(i+intervalo(1)-1)-samp(1)+1+round(maxQRSoff/1000*messages.setup.wavedet.freq))),scale)];
                end
                if   size(pontosnew,1)>1
                    weight=ones(length(pontosnew),1);
                    
                    if i==1
                        if size(timenew,2)>1 %4JUL08
                            intervalonewlead=1:timenew(i+1);
                        else
                            intervalonewlead=1:length(w1);
                        end
                    elseif i<length(timenew)
                        intervalonewlead=timenew(i-1):timenew(i+1);
                    else
                        intervalonewlead=timenew(i-1):length(w1);
                    end
                    [a2,v2] = optimline(pontosnew, weight,OPT); %#ok<ASGLU>
                else
                    v2=NaN;
                end
                if isnan(v2)
                    picoffall(end)=NaN;
                    recursioncount.QRSoff((i+intervalo(1)-1))=recursioncount.QRSoff((i+intervalo(1)-1))-1;
                else
                    if size(w1,2)==3
                        newleadbeat12=(w1(intervalonewlead,1).*v2(1)+w2(intervalonewlead,1)*v2(2)+w3(intervalonewlead,1)*v2(3))./norm(v2);
                        newleadbeat22=(w1(intervalonewlead,2).*v2(1)+w2(intervalonewlead,2)*v2(2)+w3(intervalonewlead,2)*v2(3))./norm(v2);
                        % signew2=(sig(intervalonewlead,1).*v2(1)+sig(intervalonewlead,2)*v2(2)+sig(intervalonewlead,3)*v2(3))./norm(v2);
                    else
                        newleadbeat12=(w1(intervalonewlead,1).*v2(1)+w2(intervalonewlead,1)*v2(2))./norm(v2);
                        newleadbeat22=(w1(intervalonewlead,2).*v2(1)+w2(intervalonewlead,2)*v2(2))./norm(v2);
                        % signew2=(sig(intervalonewlead,1).*v2(1)+sig(intervalonewlead,2)*v2(2))./norm(v2);
                    end
                    waux(intervalonewlead,1)=newleadbeat12;
                    waux(intervalonewlead,2)=newleadbeat22;
                    positionaux=position;
                    positionaux.R(i+intervalo(1)-1)=position_safe.R(i+intervalo(1)-1); %13.AGO.07
                    [positionaux,qrspiconset2,qrspicoffset2]= qrswavef(heasig,samp,time0(i),positionaux,waux,[i+intervalo(1)-1 i+intervalo(1)-1],messages); %#ok<ASGLU>
                    picoffall=[picoffall qrspicoffset2]; %#ok<AGROW>
                    amppicoffall=[amppicoffall (newleadbeat22(qrspicoffset2-intervalonewlead(1)+1))]; %#ok<AGROW>
                    amppicoffaux= newleadbeat22aux(qrspicoffset2-intervalonewlead(1)+1);
                    newleadbeat22aux=newleadbeat22;
                    newQRSoff=[newQRSoff position.QRSoff(i+intervalo(1)-1)]; %#ok<AGROW>
                    % stop criteria iii) | iv) => iv)  ou novo
                    %if  (positionaux.QRSoff(i+intervalo(1)-1)<=position.R(i+intervalo(1)-1) &...
                    %position.QRSoff(i+intervalo(1)-1)>=position.R(i+intervalo(1)-1) ) |((abs(amppicoffall(end))<abs(amppicoffall(end-1)))&...
                    %((amppicoffaux(end)*amppicoffall(end))>0) & ~isnan(position.S(i+intervalo(1)-1)))|(isnan(picoffall(end)))
                    %13.AGO.07
                    
                    % stop criteria iii) | iv) => iv) %opçao 0 parte da alternativa1e2
                    if ((abs(amppicoffall(end))<abs(amppicoffall(end-1)))&&((amppicoffaux(end)*amppicoffall(end))>0) &&...
                            ~isnan(position.S(i+intervalo(1)-1)))||(isnan(picoffall(end)))
                        %alternative criteria stop criteria iii) without S condition | iv) => iv) alternativa3
                        %if ((abs(amppicoffall(end))<abs(amppicoffall(end-1)))&((amppicoffaux(end)*amppicoffall(end))>0) )|(isnan(picoffall(end)))
                        %alternative criteria stop criteria  iv) => iv)
                        %if (isnan(picoffall(end)))
                        if isnan(picoffall(end))
                            position.QRSoffcriteria(i+intervalo(1)-1)=4;
                        else
                            position.QRSoffcriteria(i+intervalo(1)-1)=3;
                        end
                        picoffall(end)=NaN;
                        recursioncount.QRSoff((i+intervalo(1)-1))=recursioncount.QRSoff((i+intervalo(1)-1))-1;
                    else
                        %if ~(abs(newQRSoff(end)-newQRSoff(end-1))>1)
                        %        position.QRSoffcriteria(i+intervalo(1)-1)=1;
                        %elseif nansum(newQRSoff(end)-newQRSoff(1:end-1)==0)>1
                        if nansum(newQRSoff(end)-newQRSoff(1:end-1)==0)>1 % alternativa2
                            position.QRSoffcriteria(i+intervalo(1)-1)=2;
                        end
                        position.QRSoff(i+intervalo(1)-1) =  positionaux.QRSoff(i+intervalo(1)-1) ;
                        position.S(i+intervalo(1)-1) = positionaux.S(i+intervalo(1)-1);
                        %position.R_inQRSon(i+intervalo(1)-1) = positionaux.R(i+intervalo(1)-1) ;
                        position_safe.R(i+intervalo(1)-1) = positionaux.R(i+intervalo(1)-1) ;
                        position.Rprima(i+intervalo(1)-1) = positionaux.Rprima(i+intervalo(1)-1);
                    end
                end
                if (abs(newQRSoff(end)-newQRSoff(end-1))>1)
                    position.QRSoffcriteria(i+intervalo(1)-1)=1;
                end
            end
        end
        
    end
    %%% protection: if 2D is being used and QRSon was marked after R
    %%% %05Sep07
    if size(w1,2)==2 && position.QRSoff(i+intervalo(1)-1)<=(1+position.qrs(i+intervalo(1)-1)-1)
        position.QRSoff(i+intervalo(1)-1)=NaN;
        position.QRSoffcriteria(i+intervalo(1)-1)=-1;
    end
end