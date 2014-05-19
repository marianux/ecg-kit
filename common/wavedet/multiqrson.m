function [position0,position,messages]= multiqrson(heasig,w1,w2,w3,auxqrspiconset,QRSon,time,timenew,position0,position,intervalo,samp,messages)
% [position0,position]= multiqrson(heasig,w1,w2,w3,auxqrspiconset,QRSon,time,timenew,position0,position,intervalo,position1,position2,position3,samp)
%
% INPUT:
%   heasig: header information
%   w1,w2,w3: WT scales 1 to 3
%   auxqrspiconset: earliest relevant modulus maximum location associated to QRS complex found in any lead
%   QRSon:earliest QRS onset location in any lead (initial QRS onset location)
%   time: beats processed in this segment in each lead
%   timenew: beats processed in this segment (multilead)
%   position: struct vector with the detected points at the first multilead
%   step
%   position: struct vector with the detected points
%   intervalo: numeration of the beats processed in this segment
%   position1,position2,position3: struct vectors with the detected points using each specified lead
%   samp: samples included in the current excerpt (borders excluded)
%
%
% OUTPUT:
% actualized position0,position
%
% Rute Almeida  14.APR.2005 (based  on multiqrs)
% Last update:  Rute %05AGO2011 

global recursioncount OPT
spf=heasig.spf_ecg; %24.03.09

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
if ~isfield(messages.setup.wavedet,'maxQRSon')
    messages.setup.wavedet.maxQRSon =6.5*2; % based on CSE tolerance
end
% if ~isfield(messages.setup.wavedet,'tolsamepeak')
%         messages.setup.wavedet.tolsamepeak=5; %%%% number in samples to consider that a different peak was chosen      
% end
%tolsamepeak=messages.setup.wavedet.tolsamepeak;

% Kron=messages.setup.wavedet.Kron;
% Kroff=messages.setup.wavedet.Kroff;
% Kq=messages.setup.wavedet.Kq;
% Ks=messages.setup.wavedet.Ks;
maxQRSon=messages.setup.wavedet.maxQRSon;

auxqrson=nanmin(QRSon-samp(1)+1)-round(maxQRSon/1000*heasig.freq*spf);
%auxqrson=min(QRSon-samp(1)+1)-round(maxQRSon/1000*heasig.freq*spf); 28Abril06
for i=1:size(time,2)
       
    position.QRSonsetcriteria((i+intervalo(1)-1))=2;% default % 1 not considered
    scale=2;
    if size(w1,2)==3
        pontos=[w1(auxqrson(i):auxqrspiconset(i),scale) w2(auxqrson(i):auxqrspiconset(i),scale) w3(auxqrson(i):auxqrspiconset(i),scale)];
    else
        pontos=[w1(auxqrson(i):auxqrspiconset(i),scale) w2(auxqrson(i):auxqrspiconset(i),scale)];
    end
    %weight=ones(length(pontos),1); %RUTE 24MAR
    weight=ones(size(pontos,1),1);
    [a,v] = optimline(pontos, weight,OPT); %#ok<ASGLU>
    
    %time0=max(time);% Rute 18.07.2007
    %time0=round(nanmean(time));% Rute 18.07.2007
    time0=round(nanmedian(time));% Rute 01.03.2010
    %waux=w1; %Rute6NOV08
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

    [position0,qrspiconset0,qrspicoffset0]= qrswavef(heasig,samp,time0(i),position0,waux,[i+intervalo(1)-1 i+intervalo(1)-1],messages); %#ok<NASGU>

    if isempty(qrspiconset0)||isnan(qrspiconset0) % RUTE 6ABRI09
%         piconall= NaN;
%         amppiconall= NaN;
%         amppiconaux= NaN;
        recursioncount.QRSon((i+intervalo(1)-1))=-2;
        position.QRSon(i+intervalo(1)-1) = NaN;
        position.Q(i+intervalo(1)-1) = NaN;
        position_safe.R(i+intervalo(1)-1) = NaN; %13.AGO.07
        position.R_inQRSon(i+intervalo(1)-1) = NaN;
    else
        if size(w1,2)==3
            pontos2=[w1((position0.QRSon(i+intervalo(1)-1)-samp(1)+1-round(maxQRSon/1000*heasig.freq*spf)):qrspiconset0,scale) w2((position0.QRSon(i+intervalo(1)-1)-samp(1)+1-round(maxQRSon/1000*heasig.freq*spf)):qrspiconset0,scale) w3((position0.QRSon(i+intervalo(1)-1)-samp(1)+1)-round(maxQRSon/1000*heasig.freq*spf):qrspiconset0,scale)];
        else
            pontos2=[w1((position0.QRSon(i+intervalo(1)-1)-samp(1)+1-round(maxQRSon/1000*heasig.freq*spf)):qrspiconset0,scale) w2((position0.QRSon(i+intervalo(1)-1)-samp(1)+1-round(maxQRSon/1000*heasig.freq*spf)):qrspiconset0,scale)];
        end
        weight=ones(length(pontos2),1);
        [a2,v2] = optimline(pontos2, weight,OPT); %#ok<ASGLU>
       
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
        %waux=w1; %Rute6NOV08
        waux=nan*ones(size(w1)); %Rute6NOV08

        waux(intervalonewlead,1)=newleadbeat12;
        waux(intervalonewlead,2)=newleadbeat22;
%         positionaux=position; %RUTE 11.07.2007
        %[position,qrspiconset2,qrspicoffset2]= qrswavef(heasig,samp,time0(i),position,waux,[i+intervalo(1)-1 i+intervalo(1)-1],Kron,Kroff,Kq,Ks);
        [positionaux,qrspiconset2,qrspicoffset2]= qrswavef(heasig,samp,time0(i),position,waux,[i+intervalo(1)-1 i+intervalo(1)-1],messages); %#ok<NASGU>

        %%%%%%%%%%%%%%%%%%%%%%%RUTE 11.07.2007
        position.QRSon(i+intervalo(1)-1) = positionaux.QRSon(i+intervalo(1)-1);
        position.Q(i+intervalo(1)-1) = positionaux.Q(i+intervalo(1)-1);
        position_safe.R(i+intervalo(1)-1) = positionaux.R(i+intervalo(1)-1) ;%13.AGO.07
        position.R_inQRSon(i+intervalo(1)-1) = positionaux.R(i+intervalo(1)-1) ;
        %%%%%%%%%%%%%%%%%%%%%%%
        %if isempty(qrspiconset2) % RUTE 3ABRI09
        if isempty(qrspiconset2)||isnan(qrspiconset2) % RUTE 3ABRI09
            piconall=[qrspiconset0 NaN];
            amppiconall=[(newleadbeat2(qrspiconset0-intervalonewlead(1)+1)) NaN];
            amppiconaux= NaN;
        else
            piconall=[qrspiconset0 qrspiconset2];
            amppiconall=[(newleadbeat2(qrspiconset0-intervalonewlead(1)+1)) (newleadbeat22(qrspiconset2-intervalonewlead(1)+1))];
            amppiconaux= newleadbeat2(qrspiconset2-intervalonewlead(1)+1);
        end
        newQRSon=[NaN position0.QRSon(i+intervalo(1)-1) position.QRSon(i+intervalo(1)-1)];
        newleadbeat22aux=newleadbeat22;
        % stop criteria iii)  | iv)) ou novo %13.AGO.07
        %if position.QRSon(i+intervalo(1)-1)>=position.R(i+intervalo(1)-1)&((position0.QRSon(i+intervalo(1)-1)<position.R(i+intervalo(1)-1)))|((abs(amppiconall(end))<abs(amppiconall(end-1)))&((amppiconaux(end)*amppiconall(end))>0)  &   ~isnan(position0.Q(i+intervalo(1)-1)))|( isnan(piconall(end))) % RUTE 9.07.2007
        % stop criteria iii)  | iv)) % opçao0 % parte da alternativa2
        if ((abs(amppiconall(end))<abs(amppiconall(end-1)))&&((amppiconaux(end)*amppiconall(end))>0)  &&...
                ~isnan(position0.Q(i+intervalo(1)-1)))||( isnan(piconall(end))) % RUTE 9.07.2007
            % stop criteria iii)  | iv | i))) % alternativa1
            %if ((abs(amppiconall(end))<abs(amppiconall(end-1)))&((amppiconaux(end)*amppiconall(end))>0)  &   ~isnan(position0.Q(i+intervalo(1)-1)))|( isnan(piconall(end))) | (abs(newQRSon(end)-newQRSon(end-1))<=1) % RUTE 9.07.2007
            %alternative criteria stop criteria iii) without Q condition|iv)alternativa3
            %if ((abs(amppiconall(end))<abs(amppiconall(end-1)))&((amppiconaux(end)*amppiconall(end))>0))|( isnan(piconall(end))) % | (abs(newQRSon(end)-newQRSon(end-1))<=1) % RUTE 9.07.2007
            %alternative criteria stop criteria just iv)
            %if ( isnan(piconall(end))) % | (abs(newQRSon(end)-newQRSon(end-1))<=1) % RUTE 9.07.2007

            if isnan(piconall(end))
                position.QRSonsetcriteria((i+intervalo(1)-1))=4;
            elseif (abs(newQRSon(end)-newQRSon(end-1))<=1)% alternativa2 - exclude criteria i)
                position.QRSonsetcriteria((i+intervalo(1)-1))=1;
            else
                position.QRSonsetcriteria((i+intervalo(1)-1))=3;
            end
            recursioncount.QRSon((i+intervalo(1)-1))=-1;
            position.QRSon(i+intervalo(1)-1) = position0.QRSon(i+intervalo(1)-1);
            position.Q(i+intervalo(1)-1) = position0.Q(i+intervalo(1)-1);
            position_safe.R(i+intervalo(1)-1) = position0.R(i+intervalo(1)-1) ; %13.AGO.07
            position.R_inQRSon(i+intervalo(1)-1) = position0.R(i+intervalo(1)-1);
        else
            recursioncount.QRSon((i+intervalo(1)-1))=0;
            % stop criteria  iv)  | ii) | i) % opçao0 % parte da alternativa1
            %while ~isnan(piconall(end)) &  nansum(newQRSon(end)-newQRSon(1:end-1)==0)<2 & (abs(newQRSon(end)-newQRSon(end-1))>1 )%%%
            %alternative criteria stop criteria iv)  | ii))  % alternativa2e3
            while ~isnan(piconall(end)) &&  nansum(newQRSon(end)-newQRSon(1:end-1)==0)<2 %& (abs(newQRSon(end)-newQRSon(end-1))>1 )%%%
                recursioncount.QRSon((i+intervalo(1)-1))=recursioncount.QRSon((i+intervalo(1)-1))+1;
                if size(w1,2)==3
                    pontosnew=[w1((position.QRSon(i+intervalo(1)-1)-samp(1)+1-round(maxQRSon/1000*heasig.freq*spf)):qrspiconset2,scale) w2((position.QRSon(i+intervalo(1)-1)-samp(1)+1-round(maxQRSon/1000*heasig.freq*spf)):qrspiconset2,scale) w3((position.QRSon(i+intervalo(1)-1)-samp(1)+1)-round(maxQRSon/1000*heasig.freq*spf):qrspiconset2,scale)];
                else
                    pontosnew=[w1((position.QRSon(i+intervalo(1)-1)-samp(1)+1-round(maxQRSon/1000*heasig.freq*spf)):qrspiconset2,scale) w2((position.QRSon(i+intervalo(1)-1)-samp(1)+1-round(maxQRSon/1000*heasig.freq*spf)):qrspiconset2,scale)];
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

                    piconall(end)=NaN;
                    recursioncount.QRSon((i+intervalo(1)-1))=recursioncount.QRSon((i+intervalo(1)-1))-1;
                else
                    if size(w1,2)==3
                        newleadbeat12=(w1(intervalonewlead,1).*v2(1)+w2(intervalonewlead,1)*v2(2)+w3(intervalonewlead,1)*v2(3))./norm(v2);
                        newleadbeat22=(w1(intervalonewlead,2).*v2(1)+w2(intervalonewlead,2)*v2(2)+w3(intervalonewlead,2)*v2(3))./norm(v2);
                       %signew2=(sig(intervalonewlead,1).*v2(1)+sig(intervalonewlead,2)*v2(2)+sig(intervalonewlead,3)*v2(3))./norm(v2);
                    else
                        newleadbeat12=(w1(intervalonewlead,1).*v2(1)+w2(intervalonewlead,1)*v2(2))./norm(v2);
                        newleadbeat22=(w1(intervalonewlead,2).*v2(1)+w2(intervalonewlead,2)*v2(2))./norm(v2);
                        %signew2=(sig(intervalonewlead,1).*v2(1)+sig(intervalonewlead,2)*v2(2))./norm(v2);
                    end
                    waux(intervalonewlead,1)=newleadbeat12;
                    waux(intervalonewlead,2)=newleadbeat22;
                    positionaux=position;
                    positionaux.R(i+intervalo(1)-1)=position_safe.R(i+intervalo(1)-1); %13.AGO.07

                    [positionaux,qrspiconset2,qrspicoffset2]= qrswavef(heasig,samp,time0(i),positionaux,waux,[i+intervalo(1)-1 i+intervalo(1)-1],messages); %#ok<NASGU>
                    if isempty(qrspiconset2)||isnan(qrspiconset2) % RUTE 3ABRI09
                        piconall=[piconall NaN]; %#ok<AGROW>
                        amppiconall=[amppiconall NaN]; %#ok<AGROW>
                            amppiconaux= NaN;
                    else
                    piconall=[piconall qrspiconset2]; %#ok<AGROW>
                    amppiconall=[amppiconall (newleadbeat22(qrspiconset2-intervalonewlead(1)+1))]; %#ok<AGROW>
                    amppiconaux= newleadbeat22aux(qrspiconset2-intervalonewlead(1)+1);
                    end

                    newleadbeat22aux=newleadbeat22;
                    newQRSon=[newQRSon position.QRSon(i+intervalo(1)-1)]; %#ok<AGROW>

                    % stop criteria iii) | iv) => iv) ou novo
                    %if (positionaux.QRSon(i+intervalo(1)-1)>=position.R(i+intervalo(1)-1)& position.QRSon(i+intervalo(1)-1)<position.R(i+intervalo(1)-1))|((abs(amppiconall(end))<abs(amppiconall(end-1)))&((amppiconaux(end)*amppiconall(end))>0) & ~isnan(position.Q(i+intervalo(1)-1)))|( isnan(piconall(end)))
                    %13.AGO.07
                    % stop criteria iii) | iv) => iv) %opçao 0 % parte de alternativa1e2
                    if ((abs(amppiconall(end))<abs(amppiconall(end-1)))&&((amppiconaux(end)*amppiconall(end))>0) &&...
                            ~isnan(position.Q(i+intervalo(1)-1)))||( isnan(piconall(end)))
                        %alternative criteria stop criteria iii) without Qcondition | iv) => iv) alternativa3
                        %if ((abs(amppiconall(end))<abs(amppiconall(end-1)))&((amppiconaux(end)*amppiconall(end))>0) )|( isnan(piconall(end)))
                        %alternative criteria stop criteria just iv) => iv)
                        %if ( isnan(piconall(end)))
                        if ( isnan(piconall(end)))
                            position.QRSonsetcriteria((i+intervalo(1)-1))=4;
                        else
                            position.QRSonsetcriteria((i+intervalo(1)-1))=3;
                        end
                        piconall(end)=NaN;
                        recursioncount.QRSon((i+intervalo(1)-1))=recursioncount.QRSon((i+intervalo(1)-1))-1;
                    else
                        % if (abs(newQRSon(end)-newQRSon(end-1))<2 )
                        %     position.QRSonsetcriteria((i+intervalo(1)-1))=1;
                        %elseif nansum(newQRSon(end)-newQRSon(1:end-1)==0)>1
                        if nansum(newQRSon(end)-newQRSon(1:end-1)==0)>1 % so em alternativa2
                            position.QRSonsetcriteria((i+intervalo(1)-1))=2;
                        end
                        position.QRSon(i+intervalo(1)-1)=positionaux.QRSon(i+intervalo(1)-1) ;
                        position.Q(i+intervalo(1)-1) = positionaux.Q(i+intervalo(1)-1);
                        position_safe.R(i+intervalo(1)-1) = positionaux.R(i+intervalo(1)-1) ;
                        position.R_inQRSon(i+intervalo(1)-1) = positionaux.R(i+intervalo(1)-1) ;
                    end                  
                end
                if (abs(newQRSon(end)-newQRSon(end-1))>1 )
                    position.QRSoncriteria(i+intervalo(1)-1)=1;
                end
            end
        end

    end
    %%% protection: if 2D is being used and QRSon was marked after R
    %%% %05Sep07
    if size(w1,2)==2 && position.QRSon(i+intervalo(1)-1)>=(-1+position.qrs(i+intervalo(1)-1)) %invalid nark

        position.QRSon(i+intervalo(1)-1)=NaN;
        position.QRSoncriteria(i+intervalo(1)-1)=-1;
    end
    position.R(i+intervalo(1)-1)=position_safe.R(i+intervalo(1)-1); %13.AGO.07 
end