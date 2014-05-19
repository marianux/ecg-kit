function     [position,janelas,messages]= twavef(heasig,samp,time,position,w,intervalo,messages)
% script which identifies T wave, as well as its onset and offset
% implementation of changes to improve FN and reduce mean error utilizando a escala 5
%
%
%Input Parameters:
%   heasig: struct vector with header information
%   samp: samples included in the current excerpt (borders excluded)
%   time:  QRS times in inedexes refering the interval included in the current excerpt (borders excluded)
%   position: struct vector with the detected points
%   w: matrix with WT scales 1 to 5
%   intervalo: numeration of the beats processed in this segment
%
%Output Parameters:
%   janelas: T wave seach windows
%   actualized parameters: position
%
% Rute Almeida
% based on twave.m by Juan Pablo Martínez Cortés
% Last update: Rute Almeida  07FEB2012
%
% Designed for MATLAB Version R12; tested with MATLAB Version R13
%
%%%%%% Constants and Thresholds  !!!!!!!!!!!!!!!!!!!!!!!!!
if ~isfield(messages.setup.wavedet,'Kton')
    messages.setup.wavedet.Kton = 4;      % 2 4!
end
if ~isfield(messages.setup.wavedet,'Ktoff')
    messages.setup.wavedet.Ktoff = 2.5;% 3.5 3!
end
if ~isfield(messages.setup.wavedet,'umbraldetT')
    messages.setup.wavedet.umbraldetT = 0.25;  % We use umbraldet*sqrt(mean(w(time(i):time(i+1),4).^2))
end
if ~isfield(messages.setup.wavedet,'umbralsig')
    messages.setup.wavedet.umbralsig  =  1/8;    % To decide if there are really 2, 1 or no peak
    %threshold to decide if there is or not a significative
    % T wave.3!
end
if ~isfield(messages.setup.wavedet,'rrant_mintol')
    messages.setup.wavedet.rrant_mintol  =  0.5;    % minimum time in sec admited for current RR to be used in T wave delineation or in Exponentially averaged RR
end
if ~isfield(messages.setup.wavedet,'rrant_maxtol')
    messages.setup.wavedet.rrant_maxtol  =  1.5;    % max time in sec admited for current RR to be used in Exponentially averaged RR
end
if ~isfield(messages.setup.wavedet,'rrant_average')
    messages.setup.wavedet.rrant_average  =  [0.8 0.2];    % Exponentially averaged RR weights
end
if ~isfield(messages.setup.wavedet,'inivent_tol')
    messages.setup.wavedet.inivent_tol  =  0.1;
end
if ~isfield(messages.setup.wavedet,'inivent_tol_S')
    messages.setup.wavedet.inivent_tol_S=0.05; % sec
end
if ~isfield(messages.setup.wavedet,'finvent_tol')
    messages.setup.wavedet.finvent_tol=0.240;% sec % it would be nioce to make it depend on round(rrmed/2)
end
if ~isfield(messages.setup.wavedet,'finvent_max')
    messages.setup.wavedet.finvent_max=0.6;% sec
end
if ~isfield(messages.setup.wavedet,'scale')
    messages.setup.wavedet.scale=4;
end
if ~isfield(messages.setup.wavedet,'scale2')
    messages.setup.wavedet.scale2=5;
end
if ~isfield(messages.setup.wavedet,'scalezerocros')
    messages.setup.wavedet.scalezerocros=3;
end
if ~isfield(messages.setup.wavedet,'scalezerocros2') %for twavetask5
    messages.setup.wavedet.scalezerocros2=4;
end

if ~isfield(messages.setup.wavedet,'min_vent')
    messages.setup.wavedet.min_vent=0.1;%sec
end
if ~isfield(messages.setup.wavedet,'Tmax_Tmin_time_min')
    messages.setup.wavedet.Tmax_Tmin_time_min=0.2;%0.15;%sec %10 ENE2012: 150 ms es poco
end
if ~isfield(messages.setup.wavedet,'Tmax_Tmin_bifasic')
    messages.setup.wavedet.Tmax_Tmin_bifasic=2.5;
end
if ~isfield(messages.setup.wavedet,'T_bound_tol')
    messages.setup.wavedet.T_bound_tol=0.12;%sec
end
Kton=messages.setup.wavedet.Kton;
Ktoff=messages.setup.wavedet.Ktoff;
umbraldetT = messages.setup.wavedet.umbraldetT;
umbralsig = messages.setup.wavedet.umbralsig;
rrant_mintol=messages.setup.wavedet.rrant_mintol;
rrant_maxtol=messages.setup.wavedet.rrant_maxtol;
rrant_average=messages.setup.wavedet.rrant_average;
inivent_tol= messages.setup.wavedet.inivent_tol;
inivent_tol_S=messages.setup.wavedet.inivent_tol_S;
finvent_tol=messages.setup.wavedet.finvent_tol;
finvent_max= messages.setup.wavedet.finvent_max;
scale=messages.setup.wavedet.scale;
scale2=messages.setup.wavedet.scale2;
scalezerocros=messages.setup.wavedet.scalezerocros;
scalezerocros2=messages.setup.wavedet.scalezerocros2; %#ok<NASGU>
min_vent=messages.setup.wavedet.min_vent;%sec
Tmax_Tmin_time_min= messages.setup.wavedet.Tmax_Tmin_time_min;%sec
Tmax_Tmin_bifasic=messages.setup.wavedet.Tmax_Tmin_bifasic;% extra criteria
T_bound_tol=messages.setup.wavedet.T_bound_tol;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(heasig)
    
end
% janelas=[];
% endwinT=[];

% Initialization of auxiliary variables
if length(time)>1
    rrant = time(2)-time(1);
else
    aux=find(position.qrs~=0);
    if length(aux)>1
        rrant= diff(position.qrs(aux(end-1):aux(end))); %18MAR09
    else
        rrant=messages.setup.wavedet.freq;
    end
end

%18MAR09 % for the case in which the first RR is very wrong AGO2011
if rrant<rrant_mintol*messages.setup.wavedet.freq || rrant>rrant_maxtol*messages.setup.wavedet.freq;
    messages.warnings=[messages.warnings {['RR lower than' rrant_mintol*messages.setup.wavedet.freq ' ms or higher than' rrant_maxtol*messages.setup.wavedet.freq ' found in twavef']}];
    rrant=rrant_mintol*messages.setup.wavedet.freq;%Jul2011
end

T = []; Tprima=[]; picon=[]; picoff=[]; Ton=[]; Toff=[]; tipoT=[];
% minapos = []; minppos=[]; maxapos=[]; maxppos=[]; mina=[]; minp=[]; maxa=[];
% maxp=[];
picon_keep=[];%multileadchange
picoff_keep=[];%multileadchange
lead_keep=[];%multileadchange



%figure
%%%plot(messages.lixo)
% hold on
janelas=NaN*ones(length(time),3);
endwinT=NaN*ones(length(time),1);
for i = 1:length(time)   % For each beat
    
    
    
    if (i>1),
        if (rrant_maxtol*rrant>(time(i)-time(i-1)))&&(time(i)-time(i-1)>rrant_mintol*rrant),
            rrmed = rrant_average(1)*rrant + rrant_average(2)*(time(i)-time(i-1));% Exponentially averaged RR new value
        else
            rrmed = rrant;  % Exponentially averaged RR old value
        end
    else                  % Only for the first in each segment of the ECG
        rrmed = rrant;
    end
    rrant = rrmed;        % For next segment
    inivent = round(inivent_tol*messages.setup.wavedet.freq);   % Begining of window
    
    
    %Rute multilead 03.Dec.04
    
    %if ~isempty(pos.S(i)),              % If there is an S wave
    %         inivent = max(inivent, pos.S(i)-pos.qrs(i)+round(0.05*messages.setup.wavedet.freq));
    %     end
    
    if ~isempty(position.S(i+intervalo(1)-1));
        inivent = max(inivent, position.S(i+intervalo(1)-1)-samp(1)+1-time(i)+round(inivent_tol_S*messages.setup.wavedet.freq));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %% changed 13/06/02 Rute
    %if rrmed >= messages.setup.wavedet.freq,  %
    %if rrmed <= messages.setup.wavedet.freq,  %   21AGO09 % for rrmed < 1 seg finiven=finvent_max seg
    if rrmed <= messages.setup.wavedet.freq && i~=length(time) && (rrant_maxtol*rrant>(time(i+1)-time(i)))&&(time(i+1)-time(i)>rrant_mintol*rrant) %02FEB2011
        finvent = round(finvent_max*messages.setup.wavedet.freq);     % End of window
    else
        finvent = round(rrmed*finvent_max);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %% changed 13/06/02 Rute
    if i~=length(time),                 % For last beat in the segment
        %finvent = min(finvent,time(i+1)-time(i)-round(finvent_tol*messages.setup.wavedet.freq));
        %finvent = min(finvent,time(i+1)-time(i)-min(round((time(i+1)-time(i))*1/3),round(finvent_tol*messages.setup.wavedet.freq)));
        finvent = min(finvent,time(i+1)-time(i)-min(round((time(i+1)-time(i))*0.4),round(finvent_tol*messages.setup.wavedet.freq)));
    elseif i~=1%Rute 18MAR09
        finvent = min(finvent,time(i)-time(i-1)-round(finvent_tol*messages.setup.wavedet.freq));
    elseif length(aux)>1
        finvent = min(finvent,diff(position.qrs(aux(end-1):aux(end)))-round(finvent_tol*messages.setup.wavedet.freq));%Rute 18MAR09
    end
    
    % We work at scale 4 (in general).
    if isempty(position.QRSoff(i+intervalo(1)-1)),   % Should never happen, but...
        begwin = min(inivent + time(i),length(w)); %Rute 4Jul2011
    else
        begwin = min(max(inivent + time(i), position.QRSoff(i+intervalo(1)-1)-samp(1)+1+1),length(w));%Rute 4Jul2011
        % if QRSoff is bad anoated it can miss T wave...
    end
    
    endwin = max(1,min(finvent + time(i),length(w))); %% Rute 20/05/02
    janelas(i,:)=[i begwin endwin];
    endwinT(i)=endwin;
    %%%%plot([begwin endwin ], messages.lixo([begwin endwin]),'*g')
    
    % Positive maxima and negative minima in the window
    maxpos = begwin + modmax(w(begwin+1:endwin,scale),2,0,+1);
    minpos = begwin + modmax(w(begwin+1:endwin,scale),2,0,-1);
    [maxim ind] = max(w(maxpos,scale));   % The biggest of the positive
    maxpos = maxpos(ind);
    [minim ind] = min(w(minpos,scale));   % The biggest of the negative
    minpos = minpos(ind);
    
    if isempty(maxpos),               % If no local positive maximum
        % The maximum will be the first
        % or the last sample
        if (w(begwin,scale)>=w(endwin,scale)) && w(begwin,scale)>0,
            maxpos = begwin; maxim = w(maxpos,scale);
        elseif (w(endwin,scale)>=w(begwin,scale)) && w(endwin,scale)>0,
            maxpos = endwin; maxim = w(maxpos,scale);
        end
    end
    if isempty(minpos),               % if no local negative minimum
        % the minimum will be the first
        % or the last sample
        if (w(begwin,scale)<=w(endwin,scale)) && w(begwin,scale)<0,
            minpos = begwin; minim = w(minpos,scale);
        elseif (w(endwin,scale)<=w(begwin,scale)) && w(endwin,scale)<0,
            minpos = endwin; minim = w(minpos,scale);
        end
    end
    
    absmax = abs(maxim);
    absmin = abs(minim);
    
    if i<length(time),
        veficaz = sqrt(mean(w(time(i):time(i+1),scale).^2));
    else                        % if last beat of the segment
        veficaz = sqrt(mean(w(time(i):end,scale).^2));
    end
    
    hay_onda = ((absmax>umbraldetT*veficaz)|(absmin>umbraldetT*veficaz));
    if isempty(hay_onda)
        hay_onda = 0;
    end
    % Rute 18/06/02
    
    if endwin-begwin<(min_vent*messages.setup.wavedet.freq) || isnan(hay_onda)% se a janela tem amplitude menor do que min_vent seg enato nao existe onda T
        hay_onda =0;
    end
    
    % Is there a wave?
    if hay_onda,
        if absmax >= absmin,    % the greatest modulus maximum is the maximum
            % Now we search the two minima nearest to maxpos, one before and one after
            minapos = max(modmax(w(begwin+1:maxpos-1,scale),2,0,-1));
            minapos = begwin + minapos; % Position of the negative minimum before the maximum
            if isempty(minapos) && (maxpos ~= begwin) && (w(begwin,scale)<0),
                minapos = begwin;% If no local minimum before the maximum, take the first sample
            end
            minppos = min(modmax(w(maxpos+1:endwin,scale),2,0,-1));
            minppos = maxpos + minppos; % Position of the positive maximum after the minimum
            if isempty(minppos) && (maxpos ~= endwin) && (w(endwin,scale)<0),
                minppos = endwin;        % If no local minimum after the maximum, take the last sample
            end
            
            mina = abs(w(minapos,scale));     % Amplitude of minimum before maximum
            minp = abs(w(minppos,scale));     % Amplitude of minimum after maximum
            
            if (mina < umbralsig*absmax), % If mina is not big enough
                mina =[];                  % forget it
                
            elseif (maxpos-minapos>Tmax_Tmin_time_min*messages.setup.wavedet.freq), %or if ther are more than Tmax_Tmin_time_min sec to maxpos
                mina = [];
            end
            if (minp < umbralsig*absmax), % If minp is not big enough
                minp =[];                  % forget it
            elseif (minppos-maxpos>Tmax_Tmin_time_min*messages.setup.wavedet.freq), % and also if there are more than Tmax_Tmin_time_min sec to maxpos
                minp =[];
            end
            
            if ~isnan(mina)&~isnan(minp), %#ok<AND2>     %%% NUEVO JP
                if (mina >= minp)&&(minp < umbralsig*absmax*Tmax_Tmin_bifasic),
                    minp = [];
                elseif (minp> mina)&&(mina < umbralsig*absmax*Tmax_Tmin_bifasic),
                    mina = [];
                end
            end
            
            % Test which modulus maxima are significative and find zero crossings
            if isempty(mina),
                if isempty(minp),
                    tipoT = 2;      % only upwards T wave
                    if maxpos - minapos > 2,         % if not !!!!!!!!!!!
                        ind = zerocros(flipud(w(minapos:maxpos,scalezerocros)));   %Scale scalezerocros !!!
                        T = maxpos - ind +1;                           % Zero crossing = T wave position
                        picoff = maxpos;                               %wavelet  peak to detect offset
                    elseif isempty(minapos);  % If there were no minimum, there is no zero crossing
                        T = picant (w(begwin:maxpos,scale),maxpos);       % Take the minimum at scale scale
                        if ~isempty(T) %%%%%%%%%%%% 14/06/02 Rute
                            picoff = maxpos;
                        else
                            picoff = []; % if did not exist a peak in scale scale there is no T
                        end %%%%%%%%%%%% 14/06/02 Rute
                    else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                        messages.warnings=[messages.warnings {'unknown case: unable to classify T wave' }];
                    end
                else      % minp exists (is significative) but mina not
                    tipoT = 0;  		%normal T wave
                    if minppos -maxpos >2,       % if not!!!!!!???
                        ind = zerocros(w(maxpos:minppos,scalezerocros));  %% 07/06/02 Rute
                        if isempty(ind) %%%%%%%%%%%%%%Rute 03/09/02
                            ind = zerocros(w(maxpos:minppos,scale));
                        end %%%%%%%%%%%%%%Rute 03/09/02
                        T = maxpos + ind -1;
                        picon = maxpos;		% For determining onset and offset
                        picoff = minppos;
                    else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                        messages.warnings=[messages.warnings {'unknown case: unable to classify T wave'}];
                    end
                end
            else
                if isempty(minp),   %  mina exists (is significative) but minp not
                    tipoT = 1;  	%inverted T wave
                    if maxpos -minapos >2,  % if not !!!!!!!!!
                        ind = zerocros(w(minapos:maxpos,scalezerocros));  % wavelet zero crossing is T wave peak  %% 07/06/02 Rute
                        if isempty(ind)%%%%%%%%%%%%%%Rute 03/09/02
                            ind = zerocros(w(minapos:maxpos,scale));  % wavelet zero crossing is T wave peak  %% 03/09/02 Rute
                        end %%%%%%%%%%%%%%Rute 03/09/02
                        T = minapos + ind -1;
                        picon = minapos;
                        picoff = maxpos;
                    else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                        messages.warnings=[messages.warnings {'unknown case: unable to classify T wave'}];
                    end
                else    % both mina and minp are significative.  Biphasic wave.
                    tipoT = 5;	% biphasic neg-pos T wave
                    if maxpos - minapos > 2,    %!!!!!!!!!!!!
                        ind = zerocros(flipud(w(minapos:maxpos,scalezerocros))); %% 07/06/02 Rute
                        if isempty(ind) %%%%%%%%%%%%%%Rute 03/09/02
                            ind = zerocros(flipud(w(minapos:maxpos,scale))); %% 03/09/02 Rute
                            
                        end%%%%%%%%%%%%%%Rute 03/09/02
                        T = maxpos - ind +1;
                        picon = minapos;
                    else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                        messages.warnings=[messages.warnings {'unknown case: unable to classify T wave'}];
                    end
                    if minppos - maxpos > 2,    %!!!!!!!!!!!!
                        ind = zerocros(flipud(w(maxpos:minppos,scalezerocros))); %% 07/06/02 Rute
                        if isempty(ind)%%%%%%%%%%%%%%Rute 03/09/02
                            %ind = zerocros(flipud(w(minapos:maxpos,scale))); %% 03/09/02 Rute %18.11.05 DUVIDA
                            %DUVIDA 18Nov05
                            ind = zerocros(flipud(w(maxpos:minppos,scale))); %% 03/09/02 Rute
                        end %%%%%%%%%%%%%%Rute 03/09/02
                        %Tprima = maxpos + ind -1;
                        %DUVIDA 18Nov05
                        Tprima = minppos- ind +1;
                        picoff = minppos;
                    else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                        messages.warnings=[messages.warnings {'unknown case: unable to classify T wave'}];
                    end
                end
            end
        else        % If the greatest modulus maximum is the minimum
            % Search two maxima, one before and one after the minimum
            maxapos = max(modmax(w(begwin+1:minpos-1,scale),2,0,1));
            maxapos = begwin + maxapos;
            if isempty(maxapos) && (minpos ~= begwin) && (w(begwin,scale)>0),
                maxapos = begwin;
            end
            maxppos = min(modmax(w(minpos+1:endwin,scale),2,0,1));
            maxppos = minpos + maxppos;
            if isempty(maxppos) && (minpos ~= endwin) && (w(endwin,scale)>0),
                maxppos = endwin;
            end
            maxa = abs(w(maxapos,scale));
            maxp = abs(w(maxppos,scale)) ;                                   % See if they are significative
            if (maxa < umbralsig*absmin)
                maxa =[];
            elseif (minpos-maxapos>Tmax_Tmin_time_min*messages.setup.wavedet.freq),
                maxa = [];
            end
            if (maxp < umbralsig*absmin),
                maxp =[];
            elseif (maxppos-minpos>Tmax_Tmin_time_min*messages.setup.wavedet.freq),
                maxp = [];
            end
            if ~isnan(maxa)&~isnan(maxp), %#ok<AND2>     %%% NUEVO JP
                if (maxa >= maxp)&&(maxp < umbralsig*absmin*Tmax_Tmin_bifasic),
                    maxp = [];
                elseif (maxp> maxa)&&(maxa < umbralsig*absmin*Tmax_Tmin_bifasic),
                    maxa = [];
                end
            end
            % Test which modulus maxima are significative and find zero crossings
            if isempty(maxa),
                if isempty(maxp),
                    tipoT = 3;      % only downwards T wave
                    if minpos - maxapos > 2, %!!!!!!!!
                        ind = zerocros(flipud(w(maxapos:minpos,scalezerocros)));   %Scale scalezerocros
                        if isempty(ind) %%%%%%%%%%%%%%Rute 03/09/02
                            ind = zerocros(flipud(w(maxapos:minpos,scale)));
                        end %%%%%%%%%%%%%%Rute 03/09/02
                        T = minpos - ind +1;
                        picoff = minpos;
                    elseif isempty(maxapos);  % If there were no maximum, there is no zero crossing.
                        T = picant (w(begwin:minpos,scale),minpos);  % minimo en escala scale.
                        if ~isempty(T) %%%%%%%%%%%% 14/06/02 Rute
                            picoff = minpos;
                        else
                            picoff = []; % if did not exist a peak in scale scale there is no T
                        end %%%%%%%%%%%% 14/06/02 Rute
                    else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                        messages.warnings=[messages.warnings {'unknown case: unable to classify T wave'}];
                    end
                else      % maxp is signficative, but not maxa
                    tipoT = 1;  %inverted T wave
                    if maxppos -minpos >2,  % !!!!!!
                        ind = zerocros(w(minpos:maxppos,scalezerocros));  %% 07/06/02 Rute
                        if isempty(ind) %%%%%%%%%%%%%%Rute 03/09/02
                            ind = zerocros(w(minpos:maxppos,scale));  %% 03/09/02 Rute
                        end %%%%%%%%%%%%%%Rute 03/09/02
                        T = minpos + ind -1;
                        picon = minpos;		% For calculating onset and offset
                        picoff = maxppos;
                    else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                        messages.warnings=[messages.warnings {'unknown case: unable to classify T wave'}];
                    end
                end
            else
                if isempty(maxp),   % maxa is significative, but not maxp
                    tipoT = 0;       %normal T wave
                    if minpos -maxapos >2,
                        ind = zerocros(w(maxapos:minpos,scalezerocros)); %% 07/06/02 Rute
                        if isempty(ind) %%%%%%%%%%%%%%%%%%%% 03/09/02 Rute
                            ind = zerocros(w(maxapos:minpos,scale)); %% 03/09/02 Rute
                        end %%%%%%%%%%%%%%%%%%%% 03/09/02 Rute
                        T = maxapos + ind -1;
                        picon = maxapos;
                        picoff = minpos;
                    else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                        messages.warnings=[messages.warnings {'unknown case: unable to classify T wave'}];
                    end
                else    % both maxa and maxp are significative.  Biphasic wave.
                    tipoT = 4;	% biphasic pos-neg T wave
                    if minpos - maxapos > 2,  %!!!!!!!!!!!
                        ind = zerocros(flipud(w(maxapos:minpos,scalezerocros))); %% 07/06/02 Rute
                        if isempty(ind) %%%%%%%%%%%%%%%%%%%% 03/09/02 Rute
                            ind = zerocros(flipud(w(maxapos:minpos,scale))); %% 03/09/02 Rute
                        end %%%%%%%%%%%%%%%%%%%% 03/09/02 Rute
                        T = minpos - ind +1;
                        picon = maxapos;
                    else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                        messages.warnings=[messages.warnings {'unknown case: unable to classify T wave'}];
                    end
                    if maxppos - minpos > 2,    %!!!!!!!!!!!!
                        ind = zerocros(flipud(w(minpos:maxppos,scalezerocros))); %% 07/06/02 Rute
                        if isempty(ind) %%%%%%%%%%%%%%%%%%%% 03/09/02 Rute
                            %ind = zerocros(flipud(w(minapos:maxpos,scale))); %% 03/09/02 Rute %18.11.05 DUVIDA
                            %DUVIDA 18Nov05
                            ind = zerocros(flipud(w(minpos:maxppos,scale))); %% 03/09/02 Rute
                        end %%%%%%%%%%%%%%%%%%%% 03/09/02 Rute
                        %Tprima = minpos + ind -1; %18.11.05 DUVIDA
                        %DUVIDA 18Nov05
                        Tprima = maxppos- ind +1;
                        picoff = maxppos;
                    else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                        messages.warnings=[messages.warnings {'unknown case: unable to classify T wave'}];
                    end
                end
            end
        end
        
        % T wave onset and offset detection
        %if isempty(T), picon=[]; picoff=[]; end
        
        picon_keep=[picon_keep picon]; %#ok<AGROW> %multileadchange
        picoff_keep=[picoff_keep picoff];  %#ok<AGROW> %multileadchange
        lead_keep=[lead_keep scale];  %#ok<AGROW> %multileadchange
        
        if ~isempty(picon),
            Ton = searchon (picon, w(max(begwin,picon-round(T_bound_tol*messages.setup.wavedet.freq)):picon,scale), Kton);
        end
        if ~isempty(picoff),
            Toff=searchoff(picoff, w(picoff:min([size(w,1) picoff+round(T_bound_tol*messages.setup.wavedet.freq) ]) ,scale) , Ktoff);
            if (Toff > endwin),
                Toff = endwin;
            end
        end
    else
        %% using scale 5 07/06/02 Rute %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Neste momento usa-se a mesma janela e os mesmos limiares que para a escala scale
        
        maxpos = begwin + modmax(w(begwin+1:endwin,scale2),2,0,+1);
        minpos = begwin + modmax(w(begwin+1:endwin,scale2),2,0,-1);
        [maxim ind] = max(w(maxpos,scale2));   % The biggest of the positive
        maxpos = maxpos(ind);
        [minim ind] = min(w(minpos,scale2));   % The biggest of the negative
        minpos = minpos(ind);
        if isempty(maxpos),               % If no local positive maximum
            % The maximum will be the first
            % or the last sample
            if (w(begwin,scale2)>=w(endwin,scale2)) && w(begwin,scale2)>0,
                maxpos = begwin; maxim = w(maxpos,scale2);
            elseif (w(endwin,scale2)>=w(begwin,scale2)) && w(endwin,scale2)>0,
                maxpos = endwin; maxim = w(maxpos,scale2);
            end
        end
        if isempty(minpos),               % if no local negative minimum
            % the minimum will be the first
            % or the last sample
            if (w(begwin,scale2)<=w(endwin,scale2)) && w(begwin,scale2)<0,
                minpos = begwin; minim = w(minpos,scale2);
            elseif (w(endwin,scale2)<=w(begwin,scale2)) && w(endwin,scale2)<0,
                minpos = endwin; minim = w(minpos,scale2);
            end
        end
        
        absmax = abs(maxim);
        absmin = abs(minim);
        if i<length(time),
            veficaz = sqrt(mean(w(time(i):time(i+1),scale2).^2));
        else                        % if last beat of the segment
            veficaz = sqrt(mean(w(time(i):end,scale2).^2));
        end
        hay_ondascale2 = ((absmax>umbraldetT*veficaz)|(absmin>umbraldetT*veficaz));
        
        if endwin-begwin<(min_vent*messages.setup.wavedet.freq) % se a janela tem amplitude menor do que min_vent seg enato nao existe onda T
            hay_ondascale2 =0;
        end
        if hay_ondascale2,
            t5; % procurar na escala scale2 10/07/02 Rute
            picon_keep=[picon_keep picon];%#ok<AGROW> %multileadchange
            picoff_keep=[picoff_keep picoff];%#ok<AGROW>%multileadchange
            lead_keep=[lead_keep scale2];%#ok<AGROW>%multileadchange
        end
    end
    % Is there a wave
    %% utilizar a escala scale2 07/06/02 Rute %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Filling the structure with positions
    if isempty(Ton), Ton = NaN; end;
    if isempty(Toff), Toff = NaN;  end;
    if isempty(T), T=NaN; end;
    if isnan(T),
        picon_keep=[picon_keep NaN]; %#ok<AGROW>%multileadchange
        picoff_keep=[picoff_keep NaN]; %#ok<AGROW>%multileadchange
        lead_keep=[lead_keep NaN]; %#ok<AGROW>%multileadchange
    end
    if isempty(Tprima), Tprima=NaN; end;
    if isempty(tipoT), tipoT=NaN; end;
    pos.Tscale(i)=scale;
    pos.Ton(i) = Ton;
    pos.Toff(i) = Toff;
    pos.T(i) = T;
    pos.Tprima(i) = Tprima;
    pos.Ttipo(i) = tipoT;
    %%%plot([Toff(~isnan(Toff))], messages.lixo([Toff(~isnan(Toff))]),'*r')
    T = []; Tprima=[]; picon=[]; picoff=[]; Ton=[]; Toff=[]; tipoT=[];
%     minapos = []; minppos=[]; maxapos=[]; maxppos=[]; mina=[]; minp=[]; maxa=[];
%     maxp=[];
    
end
position.Tscale(intervalo(1):intervalo(2)) = pos.Tscale;
position.Ton(intervalo(1):intervalo(2)) = pos.Ton+samp(1)-1;
position.Toff(intervalo(1):intervalo(2))= pos.Toff+samp(1)-1;
position.T(intervalo(1):intervalo(2)) = pos.T+samp(1)-1;
position.Tprima(intervalo(1):intervalo(2)) = pos.Tprima+samp(1)-1;
position.Ttipo(intervalo(1):intervalo(2)) = pos.Ttipo;