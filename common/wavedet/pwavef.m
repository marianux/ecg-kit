function     [position,picon_all,picoff_all,messages]= pwavef(heasig,samp,time,position,w,intervalo,ultimo_anot,messages)
% script which identifies P wave, as well as its onset and offset
% admits 4 P wave morphologies: normal, inverted, biphasic +- and biphasic -+
% looking on the scale n=4 first and then if not found looking on the scale n=5
% scale m=3 for zero crossings
% if there is not P peak, there are not onset and offset either
% if there is no zero croosing in the scale m, look for it in scale n
%
%Input Parameters:
%   heasig: struct vector with header information
%   samp: samples included in the current excerpt (borders excluded)
%   time:  QRS times in inedexes refering the interval included in the current excerpt (borders excluded)
%   position: struct vector with the detected points
%   w: matrix with WT scales 1 to 5
%   intervalo: numeration of the beats processed in this segment
%   ultimo_anot:
%
%Output Parameters:
%   actualized parameters: position
%
% Rute Almeida 
% based on pwave.m by Juan Pablo Martínez Cortés
% Last update: Rute Almeida 07FEB2012
%
% Designed for MATLAB Version R12; tested with MATLAB Version R13
%
% Initialization of auxiliary variables
P = [];  picon=[]; picoff=[]; Pon=[]; Poff=[]; % alterado a 11/09/02
Pprima=[]; %Rute 28/06/02
picon_all=[];picoff_all=[];
% changes=[1 2 3 4 5]; % codes of the new strategies to apply % Rute 08/08/02
changes=[1 2 4 5]; % 3 look for peak on scale 3 instead of scale n

if nargin<8
    messages.warnings=[];
elseif nargin<7
    messages.errors=[messages.errors {'Fatal error in wavedet_3D\pwavef: not enough inputs.'}];
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
if isempty(heasig)
end

%%%%%% Constants and Thresholds  !!!!!!!!!!!!!!!!!!!!!!!!!
if ~isfield(messages.setup.wavedet,'umbraldetP')
    messages.setup.wavedet.umbraldetP = 0.02;%0.15;  % umbraldet*sqrt(mean(w(time(i):time(i+1),4).^2))
end
if ~isfield(messages.setup.wavedet,'umbralsig')
    messages.setup.wavedet.umbralsig = 1/8; %Rute 28/06/02 %threshold to decide if there is or not a significative P wave.
end
if ~isfield(messages.setup.wavedet,'inivent_tol_P')
    messages.setup.wavedet.inivent_tol_P  =  0.34;
end
if ~isfield(messages.setup.wavedet,'finvent_tol_P')
    messages.setup.wavedet.finvent_tol_P  =  0.05;  % 0.1 it seems to be larger
end
if ~isfield(messages.setup.wavedet,'Kpon')
    messages.setup.wavedet.Kpon = 2;% 3 ;         
end
if ~isfield(messages.setup.wavedet,'Kpoff')
    messages.setup.wavedet.Kpoff=  1.1;         % 1.35;
end

if ~isfield(messages.setup.wavedet,'P_bound_tol')
    messages.setup.wavedet.P_bound_tol=0.1;%sec
end
if ~isfield(messages.setup.wavedet,'min_vent_P')
    messages.setup.wavedet.min_vent_P=0.1;%sec
end

umbraldetP=messages.setup.wavedet.umbraldetP; %#ok<NASGU>  It is used in pwave_task1 script
umbralsig=messages.setup.wavedet.umbralsig;
Kpon= messages.setup.wavedet.Kpon;
Kpoff=messages.setup.wavedet.Kpoff;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for i = 1:length(time)
if ~isempty(intervalo)
    for i = 1:length(unique(intervalo(1):intervalo(end)))
        maxapos=[];
        minapos=[];
        minppos=[];
        maxppos=[];
        % First know (or guess) were is the QRS onset for this beat    

        if ~isnan(position.QRSon(i+intervalo(1)-1)) && (position.QRSon(i+intervalo(1)-1)-samp(1)+1)<time(i)
            qrson = position.QRSon(i+intervalo(1)-1)-samp(1)+1;
        else
            qrson = time(i) - round(messages.setup.wavedet.finvent_tol_P*messages.setup.wavedet.freq);
        end

        %inivent = round(0.2*messages.setup.wavedet.freq);
        inivent = round(messages.setup.wavedet.inivent_tol_P*messages.setup.wavedet.freq);
        if (i-1+intervalo(1)-1)>1,   % Protection, window to begin after Toff of last beat
            if ~isnan(position.Toff(i-1+intervalo(1)-1)),
                inivent = min(inivent, qrson-(position.Toff(i-1+intervalo(1)-1)-samp(1)+1));
            elseif ~isnan(position.T(i-1+intervalo(1)-1)),                          % Rute 28/01/03
                inivent = min(inivent, qrson-(position.T(i-1+intervalo(1)-1)-samp(1)+1));       % Rute 28/01/03
            elseif ~isnan(position.QRSoff(i-1+intervalo(1)-1)),                     % Rute 28/01/03
                inivent = min(inivent, qrson-(position.QRSoff(i-1+intervalo(1)-1)-samp(1)+1));  % Rute 28/01/03
            else                                                % Rute 28/01/03
                inivent = min(inivent, round(qrson-(position.qrs(i-1+intervalo(1)-1)-samp(1)+1)));     % Rute 8/10/08
            end
        end 
       finvent = round(messages.setup.wavedet.finvent_tol_P*messages.setup.wavedet.freq);
        ultimo_anot = ultimo_anot - samp(1)+1;
        % last anotation of the previous segment can overlap with the present one
        if (i==1 && ultimo_anot>=1),   % Protection in case it is the first beat
            begwin = round(max([1 qrson-inivent+1 ultimo_anot+1]));
        else
            begwin = round(max(1, qrson-inivent+1));
        end
        endwin =  round(max(1,qrson- finvent+1)); %27SEP2011
        if endwin-begwin>=messages.setup.wavedet.min_vent_P*messages.setup.wavedet.freq; %27SEP2011
            n=4;    % We work at scale 4 first
            pwave_task1 % subtask finding extreme points and deciding if there is wave is the scale n
            if ~hay_onda,   % look in the scale 5
                n=5;
                pwave_task1 % subtask finding extreme points and deciding if there is wave is the scale n
            end
        else
            hay_onda=0;
        end
        % deciding the morphology and marking the peak(s), onset and offset on the scale n
        if hay_onda && (sum(changes==2)==1 || n==4)

            if sum(changes==3)==1
                m=3; % zero crossings in the m scale
            else
                m=n; % zero crossings in the m=n scale
            end
            if sum(changes==1)==1%%%%%look for other significant extreme points - biphasic P waves %%%%%%%%%%%Rute 28/06/02
                if absmax>absmin % look for a minimum in opposite side to absmin
                    if maxpos < minpos, % look for a minimum before maxpos
                        minapos = begwin + modmax(w(begwin+1:maxpos,n),2,0,-1); % minimum of minima
                        [minaim ind] = min(w(minapos,n));
                        minapos = minapos(ind);
                        if abs(minaim)<umbralsig*absmax; % non significant
                            minapos=[];
                        end
                    else %look for a minimum after maxpos
                        minppos = maxpos + modmax(w(maxpos+1:endwin,n),2,0,-1); % minimum of minima
                        [minpim ind] = min(w(minppos,n));
                        minppos = minppos(ind);
                        if abs(minpim)<umbralsig*absmax % non significant
                            minppos=[];
                        end
                    end
                else % look for a maximum in opposite side to absmax
                    if minpos<maxpos, % look for a maximum before minpos
                        maxapos = begwin + modmax(w(begwin+1:minpos,n),2,0,1); % maximum of maxima
                        [maxaim ind] = max(w(maxapos,n));
                        maxapos = maxapos(ind);
                        if abs(maxaim)<umbralsig*absmin % non significant
                            maxapos=[];
                        end
                    else % look for a maximum after minpos
                        maxppos = minpos + modmax(w(minpos+1:endwin,n),2,0,1);% maximum of maxima
                        [maxpim ind] = max(w(maxppos,n));
                        maxppos = maxppos(ind);
                        if abs(maxpim)<umbralsig*absmin  % non significant
                            maxppos=[];
                        end
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Rute 28/06/02
            if isempty(maxapos)&& isempty(minapos)&&isempty(maxppos)&& isempty(minppos) %Rute 28/06/02
                if maxpos < minpos %normal P wave
                    Ptipo=1;
                    if minpos - maxpos >2,
                        ind = zerocros(w(maxpos:minpos,m));   % First zero crossing after maxpos
                        if isempty(ind) && sum(changes==5)==1 % look in  the n scale %Rute 29/07/02
                            ind = zerocros(w(maxpos:minpos,n));
                        end
                        P = maxpos + ind -1;
                        picon = maxpos;                       % For detecting onset and offset
                        picoff = minpos;
                    else %Rute 28/06/02
                       messages.warnings=[messages.warnings {'unknown case: unable to classify P wave'}];   
                    end
                else  %inverted P wave
                    Ptipo=0;
                    if maxpos -minpos >2,
                        ind = zerocros(w(minpos:maxpos,m));   % First zero crossing after minpos
                        if isempty(ind) && sum(changes==5)==1 % look in  the n scale %Rute 29/07/02
                            ind = zerocros(w(minpos:maxpos,n));
                        end
                        P = minpos + ind -1;
                        picon = minpos;                       % For detecting onset and offset
                        picoff = maxpos;
                    else %Rute 28/06/02
                       messages.warnings=[messages.warnings {'unknown case: unable to classify P wave'}];  
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%% Rute 22/07/02 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif ~isempty(maxapos)|| ~isempty(maxppos) % P bifasica + -
                Ptipo=4;
                [extra]= sort([minpos maxpos maxapos maxppos]);
                if extra(2)-extra(1) >2,
                    ind = zerocros(w(extra(1):extra(2),m));   % First zero crossing after picon
                    if isempty(ind) && sum(changes==5)==1 % look in  the n scale %Rute 29/07/02
                        ind = zerocros(w(extra(1):extra(2),n));
                    end
                    P = extra(1) + ind -1;
                    picon=extra(1); % For detecting onset, offset
                else
                       messages.warnings=[messages.warnings {'unknown case: unable to classify P wave'}];  
                end
                if extra(3)-extra(2) >2,
                    ind = zerocros(w(extra(2):extra(3),m));   % First zero crossing after extra
                    if isempty(ind) && sum(changes==5)==1 % look in  the n scale %Rute 29/07/02
                        ind = zerocros(w(extra(2):extra(3),n));
                    end
                    Pprima = extra(2) + ind -1;
                    picoff=extra(3); % For detecting onset, offset
                else %Rute 28/06/02
                       messages.warnings=[messages.warnings {'unknown case: unable to classify P wave'}];  
                end
            elseif ~isempty(minapos)|| ~isempty(minppos) % P bifasica - +
                Ptipo=5;
                [extra]= sort([minpos maxpos minapos minppos]); % For detecting onset, offset and intermediate extreme
                if extra(2)-extra(1) >2,
                    ind = zerocros(w(extra(1):extra(2),m));   % First zero crossing after picon
                    if isempty(ind) && sum(changes==5)==1 % look in  the n scale %Rute 29/07/02
                        ind = zerocros(w(extra(1):extra(2),n));
                    end
                    P = extra(1) + ind -1;
                    picon=extra(1); % For detecting onset, offset
                else %Rute 28/06/02
                       messages.warnings=[messages.warnings {'unknown case: unable to classify P wave'}];  
                end
                if extra(3)-extra(2) >2,
                    ind = zerocros(w(extra(2):extra(3),m));   % First zero crossing after extra
                    if isempty(ind) && sum(changes==5)==1 % look in  the n scale %Rute 29/07/02
                        ind = zerocros(w(extra(2):extra(3),n));
                    end
                    Pprima = extra(2) + ind -1;
                    picoff=extra(3); % For detecting onset, offset
                else %Rute 28/06/02
                       messages.warnings=[messages.warnings {'unknown case: unable to classify P wave'}];  
                end
            else
                messages.warnings=[messages.warnings {'unknown case: unable to classify P wave'}];  
                Ptipo=NaN;
            end  %%%%%%%%%%%%%%%%%%%%%%%%%%%% Rute 22/07/02 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % P wave onset and offset detection
            if ~isempty(picon),
                Pon = searchon (picon, w(max(begwin,picon-round(messages.setup.wavedet.P_bound_tol*messages.setup.wavedet.freq)):picon,n), Kpon);
            end
            if ~isempty(picoff),
                Poff=searchoff(picoff, w(picoff:min(picoff+round(messages.setup.wavedet.P_bound_tol*messages.setup.wavedet.freq),endwin),n), Kpoff);
            end
            if isempty(P)&& sum(changes==4)==1
                Pon=[];Poff=[];P=[]; % Rute 25/07/02
            end
            if (isempty(Pon) || isempty(Poff)),
                Pon=[];Poff=[];P=[]; %  Why?!!!!!!!!!
            end
            maxapos=[]; maxppos=[];minapos=[];minppos=[];  %#ok<NASGU>
            maxpos=[]; minpos=[];
        else
            maxpos=[]; minpos=[];
        end
        % If maxima too small or too separated, there is no P wave
        % Filling the structure with positions
        if isempty(Pon), Pon = NaN; end;
        if isempty(Poff), Poff = NaN; end;
        if isempty(P), Ptipo=NaN; end;
        if isempty(P), P=NaN; n=NaN; end;
        if isempty(Pprima), Pprima=NaN; end; % Rute 22/07/02
        if isempty(picoff), picoff=NaN; end; % Rute 23/02/10
        if isempty(picon), picon=NaN; end; % Rute 23/02/10
        pos.Pon(i) = Pon;
        pos.Poff(i) = Poff;
        pos.P(i) = P;
        pos.Pprima(i) = Pprima; % Rute 22/07/02
        picon_all=[picon_all picon]; %#ok<AGROW>
        picoff_all=[picoff_all picoff]; %#ok<AGROW>
        P = [];  picon=[]; picoff=[]; Pon=[]; Poff=[];
        Pprima=[];% Rute 22/07/02
        pos.Pscale(i) = n;
        pos.Ptipo(i)=Ptipo;
    end

    if exist('pos','var')
        position.Pon(intervalo(1):intervalo(2)) = pos.Pon+samp(1)-1;
        position.Poff(intervalo(1):intervalo(2))=  pos.Poff+samp(1)-1;
        position.P(intervalo(1):intervalo(2)) = pos.P+samp(1)-1;
        position.Pprima(intervalo(1):intervalo(2)) = pos.Pprima+samp(1)-1; % Rute 22/07/02
        position.Pscale(intervalo(1):intervalo(2)) = pos.Pscale; % Rute 22/07/02
        position.Ptipo(intervalo(1):intervalo(2)) = pos.Ptipo; % Rute 22/07/02
    elseif ~isempty(intervalo)
        position.Pon(intervalo(1):intervalo(2)) = NaN;
        position.Poff(intervalo(1):intervalo(2))=  NaN;
        position.P(intervalo(1):intervalo(2)) = NaN;
        position.Pprima(intervalo(1):intervalo(2)) = NaN; % Rute 22/07/02
        position.Pscale(intervalo(1):intervalo(2)) =NaN;
        position.Ptipo(intervalo(1):intervalo(2))=NaN;
    end
end