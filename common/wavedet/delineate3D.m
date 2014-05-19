function [position0,position,messages]= delineate3D(w1,w2,w3,intreg,timenew,position0,position,intervalo,position1,position2,position3,indexes,samp,sig,messages)
%
%
% [position0,position,A,V,A2,V2]= delineate3D(w1,w2,w3,intreg,scale,timenew,position0,position,intervalo,position1,samp)
% 3D multilead delineation of T wave
%
% Input Parameters:
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
% Rute Almeida: trabalho de doutoramento 2.Dec.2004
% Last update: 05AGO2011
%
% MATLAB Version 6.5.0.180913a (R13)

global recursioncount OPT
A=[];A2=[];
V=[];V2=[];
if isempty(sig)
    
end
% Aon=[];Von=[];
% Vg=[];
% if nargin==15
% messages = new_empty_mesg;
% end
if ~isfield(messages.setup.wavedet,'rrant_mintol')
    messages.setup.wavedet.rrant_mintol  =  0.5;    % minimum time in sec admited for current RR to be used in T wave delineation or in Exponentially averaged RR
end
if ~isfield(messages.setup.wavedet,'rrant_maxtol')
    messages.setup.wavedet.rrant_maxtol  =  1.5;    % max time in sec admited for current RR to be used in Exponentially averaged RR
end
if ~isfield(messages.setup.wavedet,'rrant_average')
    messages.setup.wavedet.rrant_average  =  [0.8 0.2];    % Exponentially averaged RR weights
end
if ~isfield(messages.setup.wavedet,'scale')
    messages.setup.wavedet.scale=4;
end
if ~isfield(messages.setup.wavedet,'scale2')
    messages.setup.wavedet.scale2=5;
end
if ~isfield(messages.setup.wavedet,'T_CSE_tol')
    messages.setup.wavedet.T_CSE_tol=30.6*2/1000; % based on CSE tolerance
end
if ~isfield(messages.setup.wavedet,'Tconvergence_crit')
    messages.setup.wavedet.Tconvergence_crit=1; % 1 maks are acepted as the same if they difer by Tconvergence_crit
end

rrant_mintol=messages.setup.wavedet.rrant_mintol;
rrant_maxtol=messages.setup.wavedet.rrant_maxtol;
rrant_average=messages.setup.wavedet.rrant_average;
scale=messages.setup.wavedet.scale;
scale2=messages.setup.wavedet.scale2;
Tconvergence_crit=messages.setup.wavedet.Tconvergence_crit;

intervalo=intervalo(1):intervalo(end);

scalei=scale*ones(1,size(intreg,1));
if ~isempty(w3)
    %scalei(sum([position1.Tscale(intervalo);position2.Tscale(intervalo);position3.Tscale(intervalo)]==5)>1)=5;
    %Rute 17Mar06
    aux=[(indexes(1,intervalo));(indexes(2,intervalo));(indexes(3,intervalo))];
    aux(1,~isnan(aux(1,:)))=position1.Tscale(aux(1,~isnan(aux(1,:))));
    aux(2,~isnan(aux(2,:)))=position2.Tscale(aux(2,~isnan(aux(2,:))));
    aux(3,~isnan(aux(3,:)))=position3.Tscale(aux(3,~isnan(aux(3,:))));
    scalei(sum(aux==5)>1)=5;
    %%%%%%
else
    scalei(sum([position1.Tscale(intervalo);position2.Tscale(intervalo)]==scale2)>1)=scale2;
end

for i=1:size(intreg,1)

    if isnan(intreg(i,1))
        position.Tscale(intervalo(i))=NaN;
        position.Toff(intervalo(i))=NaN;
        position.T(intervalo(i))=NaN;
        position.Ttipooff(intervalo(i))=NaN;
        position.Ttipo(intervalo(i))=NaN;
        position.Tprima(intervalo(i))=NaN;
        position.Ton(intervalo(i))=NaN;
        position.Ttipoon(intervalo(i))=NaN;
    else

        scale=scalei(i);
        %marks taken as the last in the 3 leads

        aux2=find(~isnan(indexes(2,intervalo(i))));
        aux3=find(~isnan(indexes(3,intervalo(i))));
        S=max([position1.S(indexes(~isnan(indexes(1,intervalo(i))),intervalo(i))) position2.S(indexes(2*aux2,intervalo(i))) position3.S(indexes(3*aux3,intervalo(i))) ]-samp(1)+1);
        QRSoff=max([position1.QRSoff(indexes(~isnan(indexes(1,intervalo(i))),intervalo(i))) position2.QRSoff(indexes(2*aux2,intervalo(i))) position3.QRSoff(indexes(3*aux3,intervalo(i)))]-samp(1)+1);

        if ~isempty(w3)
            pontos=[w1(intreg(i,1):min(intreg(i,2),length(w1)),scale) w2(intreg(i,1):min(intreg(i,2),length(w2)),scale) w3(intreg(i,1):min(intreg(i,2),length(w3)),scale) ];
        else
            pontos=[w1(intreg(i,1):min(intreg(i,2),length(w1)),scale) w2(intreg(i,1):min(intreg(i,2),length(w2)),scale)];
        end
        weight=ones(size(pontos,1),1);
        [a,v] = optimline(pontos, weight,OPT);

        A=[A;a]; %#ok<AGROW>
        V=[V;v]; %#ok<AGROW>
    

        %4FEB2010
        if length(intervalo)>1
            rrant=position.qrs(intervalo(2))-position.qrs(intervalo(1));
        elseif intervalo(1)>1
            rrant=position.qrs(intervalo(1))-position.qrs(intervalo(1)-1);
        else
            rrant=messages.setup.wavedet.freq/2; %JUL2011
        end

        if (i>1),
            if (rrant_maxtol*rrant>(timenew(i)-timenew(i-1)))&&(timenew(i)-timenew(i-1)>rrant_mintol*rrant),
                rrmed = rrant_average(1)*rrant + rrant_average(2)*(timenew(i)-timenew(i-1));
            else
                rrmed = rrant;  % Exponentially averaged RR
            end
        else                  % Only for the first in each segment of the ECG
            rrmed = rrant;
        end
        rrant = rrmed;        %#ok<NASGU> % For next segment

        %WT projected: from the previous QRS complex to the following QRS complex
        if ~isempty(w3)
            if i<length(timenew),
                newleadbeat=(w1(timenew(i):timenew(i+1),scale).*v(1)+w2(timenew(i):timenew(i+1),scale)*v(2)+w3(timenew(i):timenew(i+1),scale)*v(3))./norm(v);
                %signew=(sig(timenew(i):timenew(i+1),1).*v(1)+sig(timenew(i):timenew(i+1),2)*v(2)+sig(timenew(i):timenew(i+1),3)*v(3))./norm(v);
            else                        % if last beat of the segment
                newleadbeat=(w1(timenew(i):end,scale).*v(1)+w2(timenew(i):end,scale)*v(2)+w3(timenew(i):end,scale)*v(3))./norm(v);
                %signew=(sig(timenew(i):end,1).*v(1)+sig(timenew(i):end,2)*v(2)+sig(timenew(i):end,3)*v(3))./norm(v);
            end
        else
            if i<length(timenew),
                newleadbeat=(w1(timenew(i):timenew(i+1),scale).*v(1)+w2(timenew(i):timenew(i+1),scale)*v(2))./norm(v);
                %signew=(sig(timenew(i):timenew(i+1),1).*v(1)+sig(timenew(i):timenew(i+1),2)*v(2))./norm(v);
            else                        % if last beat of the segment
                newleadbeat=(w1(timenew(i):end,scale).*v(1)+w2(timenew(i):end,scale)*v(2))./norm(v);
                %signew=(sig(timenew(i):end,1).*v(1)+sig(timenew(i):end,2)*v(2))./norm(v);
            end
        end
        [pos,picon,picoff,janelas1,messages]=twave3D( newleadbeat,timenew,i,messages.setup.wavedet.freq,S,QRSoff,samp,rrmed,messages); %#ok<ASGLU>
%         janelas1=janelas1+timenew(i)-1;
        position0.Tscale(intervalo(i))=scale;
        position0.Ton(intervalo(i))=pos.Ton;
        position0.Toff(intervalo(i))=pos.Toff;
        position0.T(intervalo(i))=pos.T;
        position0.Ttipo(intervalo(i))=pos.Ttipo;
        position0.Tprima(intervalo(i))=pos.Tprima;
        
        %position0.Wavedet(intervalo(i))=newleadbeat; % derivada wavedet Maikel 30 marzo 2009
        %position0.direccion(c)=V; % vector direcci�n Maikel 30 marzo 2009
        %pontos2=[w1(picoff:intreg(i,2),scale) w2(picoff:intreg(i,2),scale) w3(picoff:intreg(i,2),scale) ];

        if ~isnan(picoff)
            if picoff>(intreg(i,2)-1)
                messages.warnings=[messages.warnings {'Tend out of the searching interval, in the fist step of multilead!'}];
                position.Tscale(intervalo(i))=NaN;
                position.Toff(intervalo(i))=NaN;
                position.T(intervalo(i))=NaN;
                position.Ttipooff(intervalo(i))=NaN;
                position.Ttipo(intervalo(i))=NaN;
                position.Tprima(intervalo(i))=NaN;
                recursioncount.Toff(intervalo(i))=NaN;
                position.contadorToff(intervalo(i))=NaN; %Rute
            else
                if (pos.Toff-samp(1)+1)==picoff
                    endaux=min([(position0.Toff(intervalo(i))-samp(1)+1+round(messages.setup.wavedet.T_CSE_tol*messages.setup.wavedet.freq)) length(w1)]); %3JUL

                else
                    endaux=min([intreg(i,2) (position0.Toff(intervalo(i))-samp(1)+1+messages.setup.wavedet.T_CSE_tol*messages.setup.wavedet.freq) ]);
                end
                if ~isempty(w3)
                    pontos2=[w1(picoff:endaux,scale) w2(picoff:endaux,scale) w3(picoff:endaux,scale)];
                    %pontos2=[w1(janelas1(2):janelas1(3),scale) w2(janelas1(2):janelas1(3),scale) w3(janelas1(2):janelas1(3),scale)];
                else
                    pontos2=[w1(picoff:endaux,scale) w2(picoff:endaux,scale)];
                end
                weight2=ones(size(pontos2,1),1);
                [a,v] = optimline(pontos2, weight2,OPT);
                A2=[A2;a]; %#ok<AGROW>
                V2=[V2;v]; %#ok<AGROW>
                if ~isempty(w3)
                    if i<length(timenew),
                        newleadbeat2=(w1(timenew(i):timenew(i+1),scale).*v(1)+w2(timenew(i):timenew(i+1),scale)*v(2)+w3(timenew(i):timenew(i+1),scale)*v(3))./norm(v);
                        %signew2=(sig(timenew(i):timenew(i+1),1).*v(1)+sig(timenew(i):timenew(i+1),2)*v(2)+sig(timenew(i):timenew(i+1),3)*v(3))./norm(v);
                    else                        % if last beat of the segment
                        newleadbeat2=(w1(timenew(i):end,scale).*v(1)+w2(timenew(i):end,scale)*v(2)+w3(timenew(i):end,scale)*v(3))./norm(v);
                        %signew2=(sig(timenew(i):end,1).*v(1)+sig(timenew(i):end,2)*v(2)+sig(timenew(i):end,3)*v(3))./norm(v);
                    end
                else
                    if i<length(timenew),
                        newleadbeat2=(w1(timenew(i):timenew(i+1),scale).*v(1)+w2(timenew(i):timenew(i+1),scale)*v(2))./norm(v);
                        %signew2=(sig(timenew(i):timenew(i+1),1).*v(1)+sig(timenew(i):timenew(i+1),2)*v(2))./norm(v);
                    else                        % if last beat of the segment
                        newleadbeat2=(w1(timenew(i):end,scale).*v(1)+w2(timenew(i):end,scale)*v(2))./norm(v);
                        %signew2=(sig(timenew(i):end,1).*v(1)+sig(timenew(i):end,2)*v(2))./norm(v);
                    end
                end
                [pos,picon2,picoff2,janelas2,messages]=twave3D(newleadbeat2,timenew,i,messages.setup.wavedet.freq,S,QRSoff,samp,rrmed,messages); %#ok<ASGLU>
%                 janelasnew=janelas2+timenew(i)-1;
                picoffall=[picoff picoff2];
                if isnan(picoff2)
                    amppicoffall=[abs(newleadbeat(picoff-timenew(i)+1)) NaN];
                else
                    amppicoffall=[abs(newleadbeat(picoff-timenew(i)+1)) abs(newleadbeat2(picoff2-timenew(i)+1))];
                end
                newToff=[position0.Toff(intervalo(i)) pos.Toff];

                if (amppicoffall(end)<amppicoffall(end-1))
                    position.contadorToff(intervalo(i))=-1;
                    %Vg(i,:)=V(end,:);
                    position.Tscale(intervalo(i))=position0.Tscale(intervalo(i));
                    position.Toff(intervalo(i))=position0.Toff(intervalo(i));
                    position.T(intervalo(i))= position0.T(intervalo(i));
                    position.Ttipooff(intervalo(i))=position0.Ttipo(intervalo(i));
                    position.Ttipo(intervalo(i))=NaN;
                    position.Tprima(intervalo(i))=position0.Tprima(intervalo(i));
                    %position.contadorToff=recursioncount.Toff;
                    position.direccionToffopt(:,intervalo(i))=V(end,:);
                    %position.Wavedet(intervalo(i))=newleadbeat2; % derivada wavedet Maikel 30 marzo 2009
                    %position.direccion(:,intervalo(i))=V2; % vector direcci�n Maikel 30 marzo 2009
                else
                    position.contadorToff(intervalo(i))=0;
                    %Vg(i,:)=V2(end,:);
                    position.direccionToffopt(:,intervalo(i))=V2(end,:);
                    %%%%%%%%%%%%%%%%%%%%%%recursion
                    while (~isnan(picoffall(end)) &&  abs(newToff(end)-newToff(end-1))>Tconvergence_crit) && (amppicoffall(end)>amppicoffall(end-1))
                        position.contadorToff(intervalo(i))=position.contadorToff(intervalo(i))+1;
                        %maxTintervalend=maxTintervalend; % teste 14
                        if (pos.Toff-samp(1)+1)==picoff
                            endaux=(pos.Toff-samp(1)+1+round(messages.setup.wavedet.T_CSE_tol*messages.setup.wavedet.freq));
                        else
                            endaux=min([intreg(i,2) (pos.Toff-samp(1)+1+round(messages.setup.wavedet.T_CSE_tol*messages.setup.wavedet.freq))]);
                        end
                        if ~isempty(w3)
                            pontosnew=[w1(picoffall(end):endaux,scale) w2(picoffall(end):endaux,scale) w3(picoffall(end):endaux,scale)];
                            %pontosnew=[w1(picoff(end):endaux,scale) w2(picoff(end):endaux,scale) w3(picoff(end):endaux,scale)];
                            %pontosnew=[w1(janelasnew(2):janelasnew(3),scale) w2(janelasnew(2):janelasnew(3),scale) w3(janelasnew(2):janelasnew(3),scale)];
                        else
                            pontosnew=[w1(picoffall(end):endaux,scale) w2(picoffall(end):endaux,scale)];
                        end

                        if   size(pontosnew,1)>1
                            weightnew=ones(size(pontosnew,1),1);
                            [a,v] = optimline(pontosnew, weightnew,OPT);
%                             if i==6
%                                 title([ecgnr '(Tend-step' num2str(2+position.contadorToff(intervalo(i)))')'])
%                             else
%                                 close
%                             end
                            A2=[A2;a]; %#ok<AGROW>
                            V2=[V2;v]; %#ok<AGROW>
                            if ~isempty(w3)
                                if i<length(timenew),
                                    newleadbeatnew=(w1(timenew(i):timenew(i+1),scale).*v(1)+w2(timenew(i):timenew(i+1),scale)*v(2)+w3(timenew(i):timenew(i+1),scale)*v(3))./norm(v);
                                    %signewnew=(sig(timenew(i):timenew(i+1),1).*v(1)+sig(timenew(i):timenew(i+1),2)*v(2)+sig(timenew(i):timenew(i+1),3)*v(3))./norm(v);
                                else                        % if last beat of the segment
                                    newleadbeatnew=(w1(timenew(i):end,scale).*v(1)+w2(timenew(i):end,scale)*v(2)+w3(timenew(i):end,scale)*v(3))./norm(v);
                                    %signewnew=(sig(timenew(i):end,1).*v(1)+sig(timenew(i):end,2)*v(2)+sig(timenew(i):end,3)*v(3))./norm(v);
                                end
                            else
                                if i<length(timenew),
                                    newleadbeatnew=(w1(timenew(i):timenew(i+1),scale).*v(1)+w2(timenew(i):timenew(i+1),scale)*v(2))./norm(v);
                                    %signewnew=(sig(timenew(i):timenew(i+1),1).*v(1)+sig(timenew(i):timenew(i+1),2)*v(2))./norm(v);
                                else                        % if last beat of the segment
                                    newleadbeatnew=(w1(timenew(i):end,scale).*v(1)+w2(timenew(i):end,scale)*v(2))./norm(v);
                                    %signewnew=(sig(timenew(i):end,1).*v(1)+sig(timenew(i):end,2)*v(2))./norm(v);
                                end
                            end

                            [pos,picon2,picoffnew,janelasnew,messages]=twave3D(newleadbeatnew,timenew,i,messages.setup.wavedet.freq,S,QRSoff,samp,rrmed,messages); %#ok<ASGLU>
%                             janelasnew=janelas2+timenew(i)-1;
                        else
                            picoffnew=NaN;
                        end
                        picoffall=[picoffall picoffnew]; %#ok<AGROW>
                        if isnan(picoffnew)
                            amppicoffall=[amppicoffall NaN]; %#ok<AGROW>
                        else
                            amppicoffall=[amppicoffall abs(newleadbeatnew(picoffnew-timenew(i)+1))]; %#ok<AGROW>
                        end
                        newToff=[newToff pos.Toff]; %#ok<AGROW>

                        if (amppicoffall(end)<amppicoffall(end-1))
                            pos.Toff=NaN; % because the lead got worse
                        end
                        if isnan(pos.Toff); %01Jun05
                            pos.Toff=newToff(end-1);
                            position.contadorToff(intervalo(i))=position.contadorToff(intervalo(i))-1;
                        else
                            %Vg(i,:)=V2(end,:);
                            position.direccionToffopt(:,intervalo(i))=V2(end,:);
                            %contador=[contador; recursioncount.Toff];
                        end
                        if  prod(picoffall(1:end-1)-picoffall(end))==0  % ciclo infinito no teste 14
                            picoffall(end)=NaN;
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    position.Tscale(intervalo(i))=scale;
                    position.Toff(intervalo(i))=pos.Toff;
                    position.T(intervalo(i))=pos.T;
                    position.Ttipooff(intervalo(i))=pos.Ttipo;
                    position.Ttipo(intervalo(i))=NaN;
                    position.Tprima(intervalo(i))=pos.Tprima;
                    %position.contadorToff=recursioncount.Toff;
                    %position.Wavedet(intervalo(i))=newleadbeatnew; % derivada wavedet Maikel 30 marzo 2009
                    %position.direccion(:,intervalo(i))=V2; % vector direcci�n Maikel 30 marzo 2009
                
                end
                %
            end
        else
            position.Tscale(intervalo(i))=NaN;
            position.Toff(intervalo(i))=NaN;
            position.T(intervalo(i))=NaN;
            position.Ttipooff(intervalo(i))=NaN;
            position.Ttipo(intervalo(i))=NaN;
            position.Tprima(intervalo(i))=NaN;
            recursioncount.Toff(intervalo(i))=NaN;
            position.contadorToff(intervalo(i))=NaN; %Rute
            %position.Wavedet(intervalo(i))=NaN;  % derivada wavedet Maikel 30 marzo 2009
            %position.direccion(:,intervalo(i))=NaN; % vector direcci�n Maikel 30 marzo 2009   
        end
        
        Tbeginauxiliar
    end
end
% position0.A(intervalo(i))=A;
% position0.A2(intervalo(i))=A2;
% position0.V(intervalo(i))=V;
% position0.V2(intervalo(i))=V2;
