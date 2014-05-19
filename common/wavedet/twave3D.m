function [pos,picon,picoff,janelas,messages]=twave3D(w,timenew,i,freq,S,QRSoff,samp,rrmed,messages)
%
%[pos,picon,picoff,janelas]=twave3D(w,timenew,i,freq,S,QRSoff,samp,rrmed)
%
% multilead T wave delineation
%
%Input Parameters:
%   w: matrix with WT scales 1 to 5
%   timenew:  QRS times in inedexes refering the interval included in the current excerpt (borders excluded)
%   i: beat number
%   freq: sampling frequency (heasig.freq)
%   S: lattest S wave position found in any lead
%   QRSoff: lattest QRS end position found in any lead
%   samp: samples included in the current excerpt (borders excluded)
%   rrmed: Exponentially averaged RR
%
%Output Parameters:
%   pos: fiducial marks structure (position)
%   picon: position of the first relevant modulos maximum in the wavelet
%   picoff: position of the last relevant modulos maximum in the wavelet
%   janelas: T wave seach windows
%
% Rute Almeida
% Last update: Rute Almeida  05AGO2011
%
% Designed for MATLAB Version R12; tested with MATLAB Version R13
%

%%%%%% Constants and Thresholds  !!!!!!!!!!!!!!!!!!!!!!!!!
if ~isfield(messages.setup.wavedet,'umbraldetT')
    messages.setup.wavedet.umbraldetT = 0.25;  % We use umbraldet*sqrt(mean(w(time(i):time(i+1),4).^2))
end
if ~isfield(messages.setup.wavedet,'umbralsig')
    messages.setup.wavedet.umbralsig  =  1/8;
end
if ~isfield(messages.setup.wavedet,'Kton')
    messages.setup.wavedet.Kton = 4;      % 2 4!
end
if ~isfield(messages.setup.wavedet,'inivent_tol')
    messages.setup.wavedet.inivent_tol  =  0.1;
end
if ~isfield(messages.setup.wavedet,'inivent_tol_S')
    messages.setup.wavedet.inivent_tol_S=0.05; % sec
end
if ~isfield(messages.setup.wavedet,'finvent_tol')
    messages.setup.wavedet.finvent_tol=0.240;% sec
end
if ~isfield(messages.setup.wavedet,'finvent_max')
    messages.setup.wavedet.finvent_max=0.6;% sec
end
if ~isfield(messages.setup.wavedet,'min_vent')
    messages.setup.wavedet.min_vent=0.1;%sec
end
if ~isfield(messages.setup.wavedet,'Tmax_Tmin_time_min')
    messages.setup.wavedet.Tmax_Tmin_time_min=0.15;%sec
end
if ~isfield(messages.setup.wavedet,'Tmax_Tmin_bifasic')
    messages.setup.wavedet.Tmax_Tmin_bifasic=2.5;
end

Kton=messages.setup.wavedet.Kton;
Ktoff=messages.setup.wavedet.Kton;
messages.warnings=[messages.warnings {'Recall that in twave3D ktoff=kton.'}];
umbraldetT = messages.setup.wavedet.umbraldetT;
umbralsig = messages.setup.wavedet.umbralsig;
inivent_tol= messages.setup.wavedet.inivent_tol;
inivent_tol_S=messages.setup.wavedet.inivent_tol_S;
finvent_tol=messages.setup.wavedet.finvent_tol;
finvent_max= messages.setup.wavedet.finvent_max;
min_vent=messages.setup.wavedet.min_vent;
Tmax_Tmin_time_min= messages.setup.wavedet.Tmax_Tmin_time_min;%sec
Tmax_Tmin_bifasic=messages.setup.wavedet.Tmax_Tmin_bifasic;% extra criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
janelas=[];%%31May05
% Initialization of auxiliary variables
T = []; Tprima=[]; picon=[]; picoff=[]; Ton=[]; Toff=[]; tipoT=[];
% minapos = []; minppos=[]; maxapos=[]; maxppos=[]; mina=[]; minp=[]; maxa=[];
% maxp=[];
picon_keep=[];%multileadchange
picoff_keep=[];%multileadchange
lead_keep=[];%multileadchange

inivent = round(inivent_tol*freq);   % Begining of window
if ~isempty(S),              % If there is an S wave
    inivent = max(inivent, S-timenew(i)+round(inivent_tol_S*freq));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %% changed 13/06/02 Rute
if rrmed <= freq,
    %if rrmed >= freq,
    finvent = round(finvent_max*freq);     % End of window
else
    finvent = round(rrmed*finvent_max);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %% changed 13/06/02 Rute

% using timenew from diffferent leads ty ca dist less tahn 0.24 sec!!!
if i~=length(timenew),                 % For last beat in the segment
    finvent = min(finvent,timenew(i+1)-timenew(i)-round(finvent_tol*freq));
else
    finvent = min(finvent,timenew(i)-timenew(i-1)-round(finvent_tol*freq));
end

% We work at scale scale (in general).
if isempty(QRSoff),   % Should never happen, but...
    % begwin = inivent + timenew(i);
    begwin = inivent;
else
    %  begwin = max(inivent + timenew(i), QRSoff+1);
    begwin = max(inivent , QRSoff+1-timenew(i));
end
%endwin = min(finvent + timenew(i),length(w)); %% Rute 20/05/02
endwin = min(finvent ,length(w)); %% Rute 20/05/02
janelas=[janelas; i begwin endwin];%31May05
if ~isempty(begwin+1:endwin ) % 22.04.05
    
    % Positive maxima and negative minima in the window
    maxpos = begwin + modmax(w(begwin+1:endwin ),2,0,+1);
    minpos = begwin + modmax(w(begwin+1:endwin ),2,0,-1);
    [maxim ind] = max(w(maxpos ));   % The biggest of the positive
    maxpos = maxpos(ind);
    [minim ind] = min(w(minpos ));   % The biggest of the negative
    minpos = minpos(ind);
    
    if isempty(maxpos),               % If no local positive maximum
        % The maximum will be the first
        % or the last sample
        if (w(begwin )>=w(endwin )) && w(begwin )>0,
            maxpos = begwin; maxim = w(maxpos );
        elseif (w(endwin )>=w(begwin )) && w(endwin )>0,
            maxpos = endwin; maxim = w(maxpos );
        end
    end
    if isempty(minpos),               % if no local negative minimum
        % the minimum will be the first
        % or the last sample
        if (w(begwin )<=w(endwin )) && w(begwin )<0,
            minpos = begwin; minim = w(minpos );
        elseif (w(endwin )<=w(begwin )) && w(endwin )<0,
            minpos = endwin; minim = w(minpos );
        end
    end
    
    absmax = abs(maxim);
    absmin = abs(minim);
    
    if i<length(timenew),
        veficaz = sqrt(mean(w .^2)); % all w 03.Dec.04 interval was restricted before
    else                        % if last beat of the segment
        veficaz = sqrt(nanmean(w.^2)); % all w 03.Dec.04 interval was restricted before
    end
    
    hay_onda = ((absmax>umbraldetT*veficaz)|(absmin>umbraldetT*veficaz));
    
    % Rute 18/06/02
    if endwin-begwin<(min_vent*freq) % se a janela tem amplitude menor do que min_vent seg enato nao existe onda T
        hay_onda =0;
    end
    
    % Is there a wave?
    if hay_onda,
        if absmax >= absmin,    % the greatest modulus maximum is the maximum
            % Now we search the two minima nearest to maxpos, one before and one after
            minapos = max(modmax(w(begwin+1:maxpos-1 ),2,0,-1));
            minapos = begwin + minapos; % Position of the negative minimum before the maximum
            if isempty(minapos) && (maxpos ~= begwin) && (w(begwin )<0),
                minapos = begwin;% If no local minimum before the maximum, take the first sample
            end
            minppos = min(modmax(w(maxpos+1:endwin ),2,0,-1));
            minppos = maxpos + minppos; % Position of the positive maximum after the minimum
            if isempty(minppos) && (maxpos ~= endwin) && (w(endwin )<0),
                minppos = endwin;        % If no local minimum after the maximum, take the last sample
            end
            
            mina = abs(w(minapos ));     % Amplitude of minimum before maximum
            minp = abs(w(minppos ));     % Amplitude of minimum after maximum
            
            if (mina < umbralsig*absmax), % If mina is not big enough
                mina =[];                  % forget it
            elseif (maxpos-minapos>Tmax_Tmin_time_min*freq), %or if ther are more than 150 ms to maxpos
                mina = [];
            end
            if (minp < umbralsig*absmax), % If minp is not big enough
                minp =[];                  % forget it
            elseif (minppos-maxpos>Tmax_Tmin_time_min*freq), % and also if there are more than 150 ms to maxpos
                minp =[];
            end
            
            if ~isnan(mina)&~isnan(minp),      %#ok<AND2> %%% NUEVO JP
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
                        ind = zerocros(flipud(w(minapos:maxpos )));   %Scale 3 !!! % also in scale 4!!!!!!! 03.Dec.04
                        T = maxpos - ind +1;                           % Zero crossing = T wave position
                        picoff = maxpos;                               %wavelet  peak to detect offset
                    elseif isempty(minapos);  % If there were no minimum, there is no zero crossing
                        T = picant (w(begwin:maxpos ),maxpos);       % Take the minimum at scale 4
                        if ~isempty(T) %%%%%%%%%%%% 14/06/02 Rute
                            picoff = maxpos;
                        else
                            picoff = []; % if did not exist a peak in scale 4 there is no T
                        end %%%%%%%%%%%% 14/06/02 Rute
                    else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                        messages.warnings=[messages.warnings {'unknown case: unable to classify T wave.'}];
                    end
                else      % minp exists (is significative) but mina not
                    tipoT = 0;  		%normal T wave
                    if minppos -maxpos >2,       % if not!!!!!!???
                        ind = zerocros(w(maxpos:minppos ));  %% 07/06/02 Rute % also in scale 4!!!!!!! 03.Dec.04
                        %                         if isempty(ind) %%%%%%%%%%%%%%Rute 03/09/02
                        %                             ind = zerocros(w(maxpos:minppos));
                        %                         end %%%%%%%%%%%%%%Rute 03/09/02
                        T = maxpos + ind -1;
                        picon = maxpos;		% For determining onset and offset
                        picoff = minppos;
                    else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                        messages.warnings=[messages.warnings {'unknown case: unable to classify T wave.'}];
                    end
                end
            else
                if isempty(minp),   %  mina exists (is significative) but minp not
                    tipoT = 1;  	%inverted T wave
                    if maxpos -minapos >2,  % if not !!!!!!!!!
                        ind = zerocros(w(minapos:maxpos));  % wavelet zero crossing is T wave peak  %% 07/06/02 Rute % also in scale 4!!!!!!! 03.Dec.04
                        %                         if isempty(ind)%%%%%%%%%%%%%%Rute 03/09/02
                        %                             ind = zerocros(w(minapos:maxpos ));  % wavelet zero crossing is T wave peak  %% 03/09/02 Rute
                        %                         end %%%%%%%%%%%%%%Rute 03/09/02
                        T = minapos + ind -1;
                        picon = minapos;
                        picoff = maxpos;
                    else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                        messages.warnings=[messages.warnings {'unknown case: unable to classify T wave.'}];
                    end
                else    % both mina and minp are significative.  Biphasic wave.
                    tipoT = 5;	% biphasic neg-pos T wave
                    if maxpos - minapos > 2,    %!!!!!!!!!!!!
                        ind = zerocros(flipud(w(minapos:maxpos ))); %% 07/06/02 Rute % also in scale 4!!!!!!! 03.Dec.04
                        %                         if isempty(ind) %%%%%%%%%%%%%%Rute 03/09/02
                        %                             ind = zerocros(flipud(w(minapos:maxpos ))); %% 03/09/02 Rute
                        %                         end%%%%%%%%%%%%%%Rute 03/09/02
                        T = maxpos - ind +1;
                        picon = minapos;
                    else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                        messages.warnings=[messages.warnings {'unknown case: unable to classify T wave.'}];
                    end
                    if minppos - maxpos > 2,    %!!!!!!!!!!!!
                        ind = zerocros(flipud(w(maxpos:minppos ))); %% 07/06/02 Rute % also in scale 4!!!!!!! 03.Dec.04
                        %                         if isempty(ind)%%%%%%%%%%%%%%Rute 03/09/02
                        %                             ind = zerocros(flipud(w(maxpos:minppos ))); %% 03/09/02 Rute
                        %                         end %%%%%%%%%%%%%%Rute 03/09/02
                        %Tprima = maxpos + ind -1; %18.11.05
                        Tprima = minppos- ind +1;
                        picoff = minppos;
                    else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                        messages.warnings=[messages.warnings {'unknown case: unable to classify T wave.'}];
                    end
                end
            end
        else        % If the greatest modulus maximum is the minimum
            % Search two maxima, one before and one after the minimum
            maxapos = max(modmax(w(begwin+1:minpos-1 ),2,0,1));
            maxapos = begwin + maxapos;
            if isempty(maxapos) && (minpos ~= begwin) && (w(begwin )>0),
                maxapos = begwin;
            end
            maxppos = min(modmax(w(minpos+1:endwin ),2,0,1));
            maxppos = minpos + maxppos;
            if isempty(maxppos) && (minpos ~= endwin) && (w(endwin )>0),
                maxppos = endwin;
            end
            maxa = abs(w(maxapos ));
            maxp = abs(w(maxppos )) ;                                   % See if they are significative
            if (maxa < umbralsig*absmin)
                maxa =[];
            elseif (minpos-maxapos>Tmax_Tmin_time_min*freq),
                maxa = [];
            end
            if (maxp < umbralsig*absmin),
                maxp =[];
            elseif (maxppos-minpos>Tmax_Tmin_time_min*freq),
                maxp = [];
            end
            if ~isnan(maxa)&~isnan(maxp),      %#ok<AND2> %%% NUEVO JP
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
                        ind = zerocros(flipud(w(maxapos:minpos )));   %Scale 3 % also in scale 4!!!!!!! 03.Dec.04
                        %                         if isempty(ind) %%%%%%%%%%%%%%Rute 03/09/02
                        %                             ind = zerocros(flipud(w(maxapos:minpos )));
                        %                         end %%%%%%%%%%%%%%Rute 03/09/02
                        T = minpos - ind +1;
                        picoff = minpos;
                    elseif isempty(maxapos);  % If there were no maximum, there is no zero crossing.
                        T = picant (w(begwin:minpos ),minpos);  % menimo en escala 4.
                        if ~isempty(T) %%%%%%%%%%%% 14/06/02 Rute
                            picoff = minpos;
                        else
                            picoff = []; % if did not exist a peak in scale 4 there is no T
                        end %%%%%%%%%%%% 14/06/02 Rute
                    else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                        messages.warnings=[messages.warnings {'unknown case: unable to classify T wave.'}];
                    end
                else      % maxp is signficative, but not maxa
                    tipoT = 1;  %inverted T wave
                    if maxppos -minpos >2,  % !!!!!!
                        ind = zerocros(w(minpos:maxppos ));  %% 07/06/02 Rute % also in scale 4!!!!!!! 03.Dec.04
                        %                         if isempty(ind) %%%%%%%%%%%%%%Rute 03/09/02
                        %                             ind = zerocros(w(minpos:maxppos ));  %% 03/09/02 Rute
                        %                         end %%%%%%%%%%%%%%Rute 03/09/02
                        T = minpos + ind -1;
                        picon = minpos;		% For calculating onset and offset
                        picoff = maxppos;
                    else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                        messages.warnings=[messages.warnings {'unknown case: unable to classify T wave.'}];
                    end
                end
            else
                if isempty(maxp),   % maxa is significative, but not maxp
                    tipoT = 0;       %normal T wave
                    if minpos -maxapos >2,
                        ind = zerocros(w(maxapos:minpos )); %% 07/06/02 Rute % also in scale 4!!!!!!! 03.Dec.04
                        %                         if isempty(ind) %%%%%%%%%%%%%%%%%%%% 03/09/02 Rute
                        %                             ind = zerocros(w(maxapos:minpos )); %% 03/09/02 Rute
                        %end %%%%%%%%%%%%%%%%%%%% 03/09/02 Rute
                        T = maxapos + ind -1;
                        picon = maxapos;
                        picoff = minpos;
                    else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                        messages.warnings=[messages.warnings {'unknown case: unable to classify T wave.'}];
                    end
                else    % both maxa and maxp are significative.  Biphasic wave.
                    tipoT = 4;	% biphasic pos-neg T wave
                    if minpos - maxapos > 2,  %!!!!!!!!!!!
                        ind = zerocros(flipud(w(maxapos:minpos ))); %% 07/06/02 Rute % also in scale 4!!!!!!! 03.Dec.04
                        %                         if isempty(ind) %%%%%%%%%%%%%%%%%%%% 03/09/02 Rute
                        %                             ind = zerocros(flipud(w(maxapos:minpos ))); %% 03/09/02 Rute
                        %                         end %%%%%%%%%%%%%%%%%%%% 03/09/02 Rute
                        T = minpos - ind +1;
                        picon = maxapos;
                    else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                        %                         [i tipoT]
                        disp('caso nao previsto; onda T nao marcada!')
                    end
                    if maxppos - minpos > 2,    %!!!!!!!!!!!!
                        ind = zerocros(flipud(w(minpos:maxppos ))); %% 07/06/02 Rute % also in scale 4!!!!!!! 03.Dec.04
                        %                         if isempty(ind) %%%%%%%%%%%%%%%%%%%% 03/09/02 Rute
                        %                             ind = zerocros(flipud(w(minapos:maxpos ))); %% 03/09/02 Rute
                        %                         end %%%%%%%%%%%%%%%%%%%% 03/09/02 Rute
                        %Tprima = minpos + ind -1; %18.11.05 DUVIDA
                        Tprima = maxppos- ind +1;
                        picoff = maxppos;
                    else %%%%%%%%%%%%%%%%%%%%% nao faz nada!!!!!!! % 12/06/02  Rute
                        messages.warnings=[messages.warnings {'unknown case: unable to classify T wave.'}];
                    end
                end
            end
        end
        
        % T wave onset and offset detection
        %if isempty(T), picon=[]; picoff=[]; end
        
        picon_keep=[picon_keep picon]; %multileadchange
        picoff_keep=[picoff_keep picoff];%multileadchange
        lead_keep=[lead_keep 4];%multileadchange
        
        if ~isempty(picon),
            Ton = searchon (picon, w(max(begwin,picon-0.12*freq):picon ), Kton);
        end
        if ~isempty(picoff),
            Toff=searchoff(picoff, w(picoff:min([size(w,1) picoff+0.12*freq ])  ) , Ktoff);
            if (Toff > endwin),
                Toff = endwin;
            end
        end
    else
        tipoT=9;
        messages.warnings=[messages.warnings {'unknown case: unable to find T wave in multilead approach using scale 4.'}];
    end       %% utilizar a escala 5 07/06/02 Rute %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
% Filling the structure with positions
if isempty(Ton), Ton = NaN; end;
if isempty(Toff), Toff = NaN;  end;
if isempty(T), T=NaN; end;
if isnan(T),
    picon_keep=[picon_keep NaN]; %multileadchange
    picoff_keep=[picoff_keep NaN];%multileadchange
    lead_keep=[lead_keep NaN];%#ok<NASGU> %multileadchange
end
if isempty(Tprima), Tprima=NaN; end;
if isempty(tipoT), tipoT=NaN; end;
pos.Ton= Ton+samp(1)-1+timenew(i)-1;
pos.Toff= Toff+samp(1)-1+timenew(i)-1;
pos.T= T+samp(1)-1+timenew(i)-1;
pos.Tprima= Tprima+samp(1)-1+timenew(i)-1;
pos.Ttipo= tipoT;
% T = []; Tprima=[]; picon=[]; picoff=[]; Ton=[]; Toff=[]; tipoT=[];
% minapos = []; minppos=[]; maxapos=[]; maxppos=[]; mina=[]; minp=[]; maxa=[];
% maxp=[];
%end
% position.Ton(intervalo(1):intervalo(2)) = pos.Ton+samp(1)-1;
% position.Toff(intervalo(1):intervalo(2))= pos.Toff+samp(1)-1;
% position.T(intervalo(1):intervalo(2)) = pos.T+samp(1)-1;
% position.Tprima(intervalo(1):intervalo(2)) = pos.Tprima+samp(1)-1;
% position.Ttipo(intervalo(1):intervalo(2)) = pos.Ttipo;
% pos.Ton
% Ton = pos.Ton(i)
% Toff= pos.Toff(i)
% T= pos.T(i)
% Tprima = pos.Tprima(i)
% Ttipo= pos.Ttipo(i);
%timenew(i)
%picon_keep
picoff=picoff_keep+timenew(i)-1;
picon=picon_keep+timenew(i)-1;
