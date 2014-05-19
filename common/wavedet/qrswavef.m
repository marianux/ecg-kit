function     [position,qrspiconset,qrspicoffset,messages]= qrswavef(heasig,samp,time,position,w,intervalo,messages)
% script which identifies q,r or s waves in each beat.
%
% Untill the end, all indexes are referred to sig and w,
% i.e. they are not global indexes but referred to the present segment.
%
%Input Parameters:
%   heasig: struct vector with header information
%   samp: samples included in the current excerpt (borders excluded)
%   time:  QRS times in inedexes refering the interval included in the current excerpt (borders excluded)
%   position: struct vector with the detected points
%   w: matrix with WT scales 1 to 5
%   intervalo: numeration of the beats processed in this segment
%   messages.setup.wavedet.Kron,Kroff,Kq,Ks,umbral,umbral1p,umbral1a,umbral2a,umbral2p,QRS_window (optional): parameters
%   messages.setup.wavedet.QRS_wtol_paa,QRS_wtol_za,QRS_wtol_pantza,QRS_wtol_ppp,QRS_wtol_zp,QRS_wtol_ppostzp (optional): parameters
%   messages.setup.wavedet.farza,farzaa,farzpp,QRSon_window,QRSoff_window, (optional): parameters
%   messages.setup.wavedet.firstwaverestriction,Qwaverestriction,Rwaverestriction,Rprimarestriction1,Rprimarestriction2,Srestriction1,Srestriction2 (optional): parameters
%  
%
%Output Parameters:
%   qrspiconset: first relevant slope associated to each QRS complex (WT extrema)
%   qrspicoffset: last relevant slope associated to each QRS complex (WT extrema)
%   actualized parameters: position
%
%Parameters details and default values:
%relative thresholds with respect to dermax to decide if other absolute maximum lines are significative: equation 3.14 R Almeida PhD thesis 
%messages.setup.wavedet.Kron = 20; 
%messages.setup.wavedet.Kroff = 14;
%messages.setup.wavedet.Kq = 15;
%messages.setup.wavedet.Ks = 8;	      % S offset (Proyeuto de Mapi) .          2.5 3 5 8;
%
%messages.setup.wavedet.umbral = 0.2; % relative threshold respecto to dermax to decide if waves posterior to qrs are significative: equation 3.13 R Almeida PhD thesis 
%messages.setup.wavedet.umbral1p = 0.09; %relative threshold respecto to dermax to decide if waves anterior to qrs are significative: equation 3.13 R Almeida PhD thesis    
%messages.setup.wavedet.umbral1a=0.06; % relative threshold respecto to dermax to decide if waves anterior to qrs are significative: equation 3.13 R Almeida PhD thesis 
%messages.setup.wavedet.umbral2a=0.15; % same as umbral1a but second peak in bifasic waves
%messages.setup.wavedet.umbral2p=0.18; % same as umbral1p but second peak in bifasic waves
%messages.setup.wavedet.QRS_window=0.1;	% half length of the initial delineation window for QRS : equation 3.12 PhD thesis
%messages.setup.wavedet.QRS_wtol_paa=0.07;% search window length for paa (peak before pa) to find main wave cracks
%messages.setup.wavedet.QRS_wtol_za=0.08;% search window length for za (zero before pa)   
%messages.setup.wavedet.QRS_wtol_pantza=0.08; search window length for za(zero before pa) and zaa (zero before pantza)
%messages.setup.wavedet.QRS_wtol_ppp=0.13;% search window length for ppp (peak after pp)
%messages.setup.wavedet.QRS_wtol_zp=0.13;% search window length for za (zero before pa) and search window length for zpp ( zero after zp)
%messages.setup.wavedet.QRS_wtol_ppostzp=0.08;%half length of the search window for ppostzp (peak after zp)   
%messages.setup.wavedet.farza=0.16;%maximum time allowed bewteen za and qrs
%messages.setup.wavedet.farzaa=0.16;%maximum time allowed bewteen zaa and za
%messages.setup.wavedet.farzp=0.16;%maximum time allowed bewteen zp and qrs  
%messages.setup.wavedet.farzpp=0.16;%maximum time allowed bewteen zp and zpp
%messages.setup.wavedet.QRSon_window=0.08; %search window length for QRS on
%messages.setup.wavedet.QRSoff_window=0.08; %search window length for QRS end
%messages.setup.wavedet.firstwaverestriction=0.12;%If qrs onset differs more than 120 ms of qrs the first not main wave is excluded and qrs onset searched again
%messages.setup.wavedet.Qwaverestriction=0.05;%If qrs onset differs more than 5 ms from Q wave, the Q wave is excluded and qrs onset searched again
%messages.setup.wavedet.Rwaverestriction=0.08;%If qrs onset differs more than 8 ms from a first wave R, the R wave is excluded and qrs onset searched again
%messages.setup.wavedet.Rprimarestriction1=0.12;%If qrs onset differs more than 120 ms of qrs the last not main wave Rprima is excluded and qrs onset searched again
%messages.setup.wavedet.Rprimarestriction2=0.05;%If qrs onset differs more than 5 ms from the last not main wave Rprima, R prima is excluded and qrs onset searched again
%messages.setup.wavedet.Srestriction1=0.27;%If qrs onset differs more than 270 ms of qrs the last not main wave S wave is excluded and qrs onset searched again
%messages.setup.wavedet.Srestriction2=0.15;%If qrs onset differs more than 150 ms from the last not main wave  S, S wave is excluded and qrs onset searched again
% 
% Last update: Rute Almeida  07FEB2012
%
% Designed for MATLAB Version R12; tested with MATLAB Version R13

if nargin<7 || ~isfield(messages,'setup') || ~isfield(messages.setup,'wavedet')% 27JUN2011
    messages.setup.wavedet=[];
end
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
if ~isfield(messages.setup.wavedet,'umbral')
    messages.setup.wavedet.umbral = 0.2;	   
% relative threshold respecto to dermax to decide if waves paa and ppp are significative 
end
if ~isfield(messages.setup.wavedet,'umbral1p')
    messages.setup.wavedet.umbral1p = 0.09;	   
% relative threshold respecto to dermax to decide if waves posterior to qrs are significative: equation 3.13 R Almeida PhD thesis 
end
if ~isfield(messages.setup.wavedet,'umbral1a')
    messages.setup.wavedet.umbral1a=0.06;
    % relative threshold respecto to dermax to decide if waves anterior to qrs are significative: equation 3.13 R Almeida PhD thesis 
end
if ~isfield(messages.setup.wavedet,'umbral2a')
    messages.setup.wavedet.umbral2a=0.15;	 % significativas.  % anterior 0.15 --0.20     
end
if ~isfield(messages.setup.wavedet,'umbral2p') % anterior 0.14 -- 0.09
    messages.setup.wavedet.umbral2p=0.18;	    
end
if ~isfield(messages.setup.wavedet,'QRS_window') % half length of the initial delineation window for QRS : equation 3.12 PhD thesis
    messages.setup.wavedet.QRS_window=0.1;	 
     % same as QRS_wtol_ppostzp_crack corresponding to
    % half length of the search window for ppostzz crack protection     
end
if ~isfield(messages.setup.wavedet,'QRS_wtol_paa') % search window length for paa (peak before pa) to find main wave cracks
    messages.setup.wavedet.QRS_wtol_paa=0.07;	    
end
if ~isfield(messages.setup.wavedet,'QRS_wtol_za') % search window length for za (zero before pa)   
    messages.setup.wavedet.QRS_wtol_za=0.08;
end
if ~isfield(messages.setup.wavedet,'QRS_wtol_pantza') %half length of the search window for pantza (peak before za)   
    messages.setup.wavedet.QRS_wtol_pantza=0.08;
    % same as QRS_wtol_pantzaa and QRS_wtol_pantza_crack corresponding to
    % half length of the search window for pantzaa and crack protection 
    % same as and QRS_wtol_zaa corresponding to search window length for zaa (zero before pantza)
end
if ~isfield(messages.setup.wavedet,'QRS_wtol_ppp') % search window length for ppp (peak after pp)
    messages.setup.wavedet.QRS_wtol_ppp=0.13;
end
if ~isfield(messages.setup.wavedet,'QRS_wtol_zp') % search window length for za (zero before pa)   
    messages.setup.wavedet.QRS_wtol_zp=0.13;
     % same as QRS_wtol_zpp corresponding to to search window length for zpp ( zero after zp)
end
if ~isfield(messages.setup.wavedet,'QRS_wtol_ppostzp') %half length of the search window for ppostzp (peak after zp)   
    messages.setup.wavedet.QRS_wtol_ppostzp=0.08;
end

if ~isfield(messages.setup.wavedet,'farza') %maximum time allowed bewteen za and qrs
    messages.setup.wavedet.farza=0.16;
end
if ~isfield(messages.setup.wavedet,'farzaa') %maximum time allowed bewteen zaa and za
    messages.setup.wavedet.farzaa=0.16;
end
if ~isfield(messages.setup.wavedet,'farzp') %maximum time allowed bewteen zp and qrs
    messages.setup.wavedet.farzp=0.16;
end
if ~isfield(messages.setup.wavedet,'farzpp') %maximum time allowed bewteen zp and qrs
    messages.setup.wavedet.farzpp=0.16;
end
if ~isfield(messages.setup.wavedet,'QRSon_window') %search window length for QRS on
    messages.setup.wavedet.QRSon_window=0.08;
end
if ~isfield(messages.setup.wavedet,'QRSoff_window') %search window length for QRS off
    messages.setup.wavedet.QRSoff_window=0.08;
end
if ~isfield(messages.setup.wavedet,'firstwaverestriction') 
    messages.setup.wavedet.firstwaverestriction=0.12;
    %If qrs onset differs more than 120 ms of qrs the first not main wave is excluded and qrs onset searched again
end
if ~isfield(messages.setup.wavedet,'Qwaverestriction') 
    messages.setup.wavedet.Qwaverestriction=0.05;
    %If qrs onset differs more than 5 ms from Q wave, the Q wave is excluded and qrs onset searched again
end
if ~isfield(messages.setup.wavedet,'Rwaverestriction') 
    messages.setup.wavedet.Rwaverestriction=0.08;
    %If qrs onset differs more than 8 ms from a first wave R, the R wave is excluded and qrs onset searched again
end
if ~isfield(messages.setup.wavedet,'Rprimarestriction1') 
    messages.setup.wavedet.Rprimarestriction1=0.12;
    %If qrs onset differs more than 120 ms of qrs the last not main wave Rprima is excluded and qrs onset searched again
end
if ~isfield(messages.setup.wavedet,'Rprimarestriction2') 
    messages.setup.wavedet.Rprimarestriction2=0.05;
%If qrs onset differs more than 5 ms from the last not main wave Rprima, R prima is excluded and qrs onset searched again

end
if ~isfield(messages.setup.wavedet,'Srestriction1') 
    messages.setup.wavedet.Srestriction1=0.27;
    %If qrs onset differs more than 270 ms of qrs the last not main wave S wave is excluded and qrs onset searched again
end
if ~isfield(messages.setup.wavedet,'Srestriction2') 
    messages.setup.wavedet.Srestriction2=0.15;
     %If qrs onset differs more than 150 ms from the last not main wave  S, S wave is excluded and qrs onset searched again
end
%%%%%%%%%%%%%%%%%%%%% Parameters !!!!!!!!!!!!!!!%%%%%%%%%%%%
QRS_window=messages.setup.wavedet.QRS_window;
QRSon_window=messages.setup.wavedet.QRSon_window; 
QRSoff_window=messages.setup.wavedet.QRSoff_window; 

QRS_wtol_paa=messages.setup.wavedet.QRS_wtol_paa;
QRS_wtol_za=messages.setup.wavedet.QRS_wtol_za;
QRS_wtol_pantza=messages.setup.wavedet.QRS_wtol_pantza;
QRS_wtol_pantza_crack=QRS_wtol_pantza;
QRS_wtol_zaa=QRS_wtol_pantza;
QRS_wtol_pantzaa=QRS_wtol_pantza;

QRS_wtol_ppp=messages.setup.wavedet.QRS_wtol_ppp;
QRS_wtol_zp=messages.setup.wavedet.QRS_wtol_zp;
QRS_wtol_ppostzp=messages.setup.wavedet.QRS_wtol_ppostzp;
QRS_wtol_ppostzp_crack=messages.setup.wavedet.QRS_window;
QRS_wtol_zpp=QRS_wtol_zp;
QRS_wtol_ppostzpp=QRS_wtol_ppostzp;

farza=messages.setup.wavedet.farza;
farzaa=messages.setup.wavedet.farzaa;
farzp=messages.setup.wavedet.farzp;
farzpp=messages.setup.wavedet.farzpp;

firstwaverestriction=messages.setup.wavedet.firstwaverestriction;	
Qwaverestriction=messages.setup.wavedet.Qwaverestriction;
Rwaverestriction=messages.setup.wavedet.Rwaverestriction;
Rprimarestriction1=messages.setup.wavedet.Rprimarestriction1;
Rprimarestriction2=messages.setup.wavedet.Rprimarestriction1;	
Srestriction1=messages.setup.wavedet.Srestriction1;	
Srestriction2=messages.setup.wavedet.Srestriction2;

umbral=messages.setup.wavedet.umbral; 
umbral1p = messages.setup.wavedet.umbral1p;    
umbral1a= messages.setup.wavedet.umbral1a;   
umbral2a = messages.setup.wavedet.umbral2a;
umbral2p = messages.setup.wavedet.umbral2p;


Kron=messages.setup.wavedet.Kron;
Kroff=messages.setup.wavedet.Kroff;
Kq=messages.setup.wavedet.Kq;
Ks=messages.setup.wavedet.Ks;
% Si a primera u a zaguera onda ye a R.  2.5 3 5 8 14
% antiguo 2.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(heasig)
    
end
% Structure for the present segment of the signal (sig) times referred to
% begining of the present segment.
pos=struct('Pon',[],'P',[],'Poff',[],'QRSon',[],'Q',[],'R',[],'Fiducial',[],'qrs',[],'Rprima',[],'S',[],'QRSoff',[],'Ton',[],'T',[],'Tprima',[],'Toff',[],'Ttipo',[],'QRSpa',[],'QRSpp',[]);

% After processing all the excerpt,pos is transferred to the struct position.
qrspiconset=[];
qrspicoffset=[];
    WT_extrema=NaN*ones(6,length(time));
    WT_zc=NaN*ones(5,length(time));
    pos_extrema=NaN*ones(6,length(time));
for i = 1:length(time), 
    qrs = time(i);
    R=[];Q=[];S=[];Rprima=[];   za=[];zaa=[];zp=[]; zpp=[]; %#ok<NASGU>
    dermax = max(abs(w(max(qrs-round(QRS_window*(messages.setup.wavedet.freq)),1):min(size(w,1),qrs+round(QRS_window*(messages.setup.wavedet.freq))),2)));
    
    pa = picant(w(max(qrs-round(QRS_window*(messages.setup.wavedet.freq)),1):qrs,2),qrs);
    % first peak before detected qrs position at scale 2
    pp = picpost(w(qrs:min(size(w,1),qrs+round(QRS_window*(messages.setup.wavedet.freq))),2),qrs);
    % first peak after detected qrs position at scale 2
    apa = w(pa,2);  % Amplitudes of pa and pp at scale 2
    app = w(pp,2);
    
    
    % If there is a crack at scale 1 but not at scale 2, it is possible that
    % pa and pp have the same sign.
    if (sign(apa)==sign(app)),        % In that case...
        paux =modmax(w(pp+1:min(pp+round(QRS_window*(messages.setup.wavedet.freq)),size(w,1)),2),2,0,-sign(app));
        if paux, pp2 = pp + paux(1);
        else pp2= []; app2=0; end %#ok<NASGU>
        % If there is the next 100 ms some modulus maximum of the contrary sign of app
        
        paux=modmax(flipud(w(max(1,pa-round(QRS_wtol_paa*(messages.setup.wavedet.freq))):pa-1,2)),2,0,-sign(apa));
        if paux, pa2 = pp - paux(1);
        else pa2=[]; end
        % If there is the next 70 ms some modulus maximum of the contrary sign of apa
        
        if isempty(pp2),
            if isempty(pa2),
                % Nothing
            else
                pa = pa2; apa = w(pa2,2);
            end
        else
            if isempty(pa2),
                pp = pp2; app = w(pp2,2);
            else
                if abs(w(pp2,2))> abs(w(pa2,2)),
                    pp = pp2;
                    app = w(pp2,2);
                else
                    pa = pa2;
                    apa = w(pa2,2);
                end
            end
        end
        
        % So, we take the peak with greater absolute value !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % Not coherent with the PFC!!!!!!!!!!!!!!!!!!
        % Think in which conditions the condition holds
        % and what happens if there is no peak ??  Solucionado esto último
    end
    
    % Next peaks before pa and after pp
    paa = modmax(flipud(w(max(1,pa-round(QRS_wtol_paa*(messages.setup.wavedet.freq))):pa-1,2)),2,umbral*dermax,sign(w(pa,2)));
    ppp = modmax(w(pp+1:fix(min(pp+round(QRS_wtol_ppp*(messages.setup.wavedet.freq)),size(w,1))),2),2,umbral*dermax,sign(w(pp,2)));
    
    % Protection against R-waves with cracks
    
    if ~isempty(paa),
        paa = pa - paa(end);
        pa = paa; 
        apa = max(w(paa,2),apa);
    end
    
    if ~isempty(ppp),
        ppp = pp + ppp(end);
        pp = ppp; 
        app = max(app,w(ppp,2));
    end
    
    % depending of the signs of apa and app we can know if the detected qrs position
    % corresponded to a positive or a negative wave.
    
    if (apa>0 & app<0),   %#ok<AND2> % qrs is positive qRs, Rsr, rsR.
        ind = zerocros(flipud(w(max(1,pa-round(QRS_wtol_za*(messages.setup.wavedet.freq))):pa-1,1)));
        za = pa -ind;         % previous zero
        if ~isempty(za),
            % search for the first negative peak before za (and thus, before the positive peak pa)
            % paux= modmax(flipud(w(za-0.08*(messages.setup.wavedet.freq):za,2)),2,0,-1);
            %  %%RUTE: can give error in i=1 if za < 0.08*(messages.setup.wavedet.freq)  %%24Mar06
            paux= modmax(flipud(w(max(1,za-round(QRS_wtol_pantza*(messages.setup.wavedet.freq))):za,2)),2,0,-1);
            if paux,  pantza = za - paux(1)+1;
            else  pantza = []; end
            apantza = w(pantza,2);
            % p ant za = peak anterior to za ;  % a p ant za = amplitude of the peak anterior to za
            
            %  Protection against waves with cracks
            %  If there is another peak of the same sign as pantza just before it, take it as the position of
            %  the peak, and as amplitude, the maximum of the two.
            %paux =modmax(flipud(w(pantza-round(0.05*(messages.setup.wavedet.freq)):pantza-1,2)),2,umbral2a*dermax,sign(w(pantza,2)));
            %  %%RUTE: can give error in i=1 if pantza < 0.08*(messages.setup.wavedet.freq)  %%24Mar06
            %paux
            %=modmax(flipud(w(max(1,pantza-round(0.05*(messages.setup.wavedet.freq))):pantza-1,2)),2,umbral2a*dermax,sign(w(pantza,2))); 
            % parameter  changed form 0.08 to 0.05 ?!?!? when? Why
            paux =modmax(flipud(w(max(1,pantza-round(QRS_wtol_pantza_crack*(messages.setup.wavedet.freq))):pantza-1,2)),2,umbral2a*dermax,sign(w(pantza,2)));
            if ~isempty(paux)
                paux = pantza -paux(end);
                pantza = paux;
                apantza = sign(w(paux,2))*max(abs(apantza), abs(w(paux,2)));
            end
            
            if (apantza>0)|(abs(apantza)<umbral1a*dermax) %#ok<OR2>   % If apantza is a little peak, forget it and za
                za = [];
            elseif (isempty(pantza))
                za =[];
            end
            
            if ~isempty(za), 	                              % If we are interested in the zero crossing "za"
                % Search for the zero crossign before "za", called "zaa"
                %ind = zerocros(flipud(w(pantza-0.08*(messages.setup.wavedet.freq):pantza-1,1)));
                %%RUTE: can give error in i=1 if pantza < 0.08*(messages.setup.wavedet.freq)
                %%14Mar06
                ind = zerocros(flipud(w(max(1,pantza-round(QRS_wtol_zaa*(messages.setup.wavedet.freq))):pantza-1,1)));  %%RUTE:  %14Mar06
                zaa = pantza-ind;    
               	
                
                if ~isempty(zaa),
                    paux= modmax(flipud(w(max(1,zaa-round(QRS_wtol_pantzaa*(messages.setup.wavedet.freq))):zaa,2)),2,0,+1); % Search positive peak before "zaa"
                    if paux,  pantzaa = zaa - paux(1)+1;
                    else  pantzaa = []; end
                    apantzaa = w(pantzaa,2);                     % Amplitude of the peak anterior to zaa: apantzaa
                    if (apantzaa<0)|(abs(apantzaa)<umbral2a*dermax), %#ok<OR2> % If little, forget it and zaa
                        zaa=[];
                    elseif(isempty(pantzaa))
                        zaa=[];
                    end
                end
            end
        end
        
        %%%%%%%%%%%%% now the same, but after pp (pico posterior) %%%%%%%%%%%%%%%%%%
        ind = zerocros(w(pp+1:min(size(w,1),pp+round(QRS_wtol_zp*(messages.setup.wavedet.freq))),1));	% Search zero after pp (zp)
        zp = pp + ind;
        if ~isempty(zp),
            paux= modmax(w(zp:min(size(w,1),zp+round(QRS_wtol_ppostzp*(messages.setup.wavedet.freq))),2),2,0,+1); % Search peak after zp (ppostzp)
            if paux,  ppostzp = zp + paux(1)-1;
            else  ppostzp = []; end
            appostzp = w(ppostzp,2);                                         % Amplitude of ppostzp at scale 2
            
            
            %  Protection against waves with cracks
            %  If there is another peak of the same sign as ppostzp just after it, take it as the position of
            %  the peak, and as amplitude, the maximum of the two.
            paux =modmax(w(ppostzp+1:min(size(w,1),ppostzp+round(QRS_wtol_ppostzp_crack*(messages.setup.wavedet.freq))),2),2,umbral2p*dermax,sign(w(ppostzp,2)));
            if ~isempty(paux)
                paux = ppostzp + paux(end);
                ppostzp = paux;
                appostzp = sign(w(paux,2))* max(abs(appostzp), abs(w(paux,2)));
            end
            
            if (appostzp<0)|(abs(appostzp)<umbral1p*dermax),     %#ok<OR2>    % if peak is not large enough
                zp=[];                                                % forget it and zp
            elseif isempty(ppostzp),
                zp=[];
            end
            if ~isempty(zp),                                         % else, search zero crossing after it
                ind = zerocros(w(ppostzp+1:min(size(w,1),ppostzp+round(QRS_wtol_zpp*(messages.setup.wavedet.freq))),1));
                zpp = ppostzp + ind;
                if ~isempty(zpp),                                      % and next peak after zero crossing
                    paux= modmax(w(zpp:min(size(w,1),zpp+round(QRS_wtol_ppostzpp*(messages.setup.wavedet.freq))),2),2,0,-1);
                    if paux,  ppostzpp = zpp + paux(1)-1;
                    else  ppostzpp = []; end
                    
                    appostzpp = w(ppostzpp,2);
                    if (appostzpp>0)|(abs(appostzpp)<umbral2p*dermax),  %#ok<OR2> % If not large enough, forget it
                        zpp=[];
                    elseif(isempty(ppostzpp))
                        zpp=[];
                    end
                end
            end
        end
        
        
        if (isempty(zaa))& (qrs-za)>round(farza*(messages.setup.wavedet.freq)), %#ok<AND2> % farza changed from 0.08 MAY2011 
            za = [];   % If first zero in the wavelet before qrs is too far from QRS
        end            % forget it
        
        if ~isempty(zaa) & (za - zaa)>round( farzaa*(messages.setup.wavedet.freq)),%#ok<AND2> % farza changed from 0.16 MAY2011 
            zaa=[];     % If zaa is too far from za, forget it
        end
        
        if isempty(zpp) & (zp-qrs)>round(farzp*(messages.setup.wavedet.freq)),%#ok<AND2> % farza changed from 0.16 MAY2011 
            zp = [];   % If zp is too far from qrs, forget it
        end
        if ~isempty(zpp) & (zpp-zp)>round(farzpp*(messages.setup.wavedet.freq)) %#ok<AND2> % farza changed from 0.16 MAY2011 
            zpp = [];  % if zpp is too far from zpp forget it
        end
        
        % wave label asignation
        if isempty(za),            % No waves before qrs (which is positive)
            R=qrs; S=zp; Rprima=zpp;  % RSR' or RS or R
            %2AGO2011
        elseif isempty(zp),        % if some wave before qrs but not waves after
            if ~isempty(zaa),         % if there are two waves befor qrs
                R=zaa; S=za; Rprima = qrs; % RSR'
            else                      % if only one wave before qrs
                Q=za; R=qrs;            % QR
            end
        else                       % if there are waves before and after qrs
            if (~isempty(zpp))&&(~isempty(zaa))   % if there are two waves before and two waves after
                if (abs(apantzaa)>abs(appostzpp)),% Forget the smaller of the two extreme waves
                    zpp=[]; 
                else zaa=[];
                end
            end
            % Now there are no more than 4 possible waves
            if ~isempty(zpp)&&(isempty(zaa)),     % if we have one wave before and two after
                if abs(appostzpp)>abs(apantza),     % if the second after is bigger than the one before
                    R = qrs; S = zp; Rprima = zpp; %	RSR'
                else                                % if the one before is bigger than the second after
                    Q = za; R = qrs; S=zp;	      % QRS
                end
            elseif ~isempty(zaa)&& isempty(zpp),   % if there are two waves before and one after
                if abs(apantzaa)>abs(appostzp),    % ...
                    R=zaa; S=za; Rprima=qrs;       %  RSR'
                else
                    Q=za; R=qrs; S=zp;             %  QRS
                end
            elseif za~=qrs && zp~=qrs                            % else = if there are one wave before and one after qrs
                Q=za; R=qrs; S=zp;             % QRS
            else %MAR06;2012
                if za==qrs %%%%% strange case R==Q???
                  if apantza<0
                     Q=za; R=zp;
                  else
                     R=za; S=zp;
                  end
                else
                    if appostzp<0
                        R=qrs; S=zp;
                    else
                        Q=qrs; R=zp;
                    end
                end
            end
            
        end
        % We begin again in case the qrs peak is negative: the process is similar
    else    % Negative QRS: Qrs qrS QS
        
        %
        ind = zerocros(flipud(w(max(1,pa-round(QRS_wtol_za*(messages.setup.wavedet.freq))):pa-1,1)));
        za = pa -ind;         % previous zero. za: zero anterior
        if ~isempty(za),
            %paux= modmax(flipud(w(za-0.08*(messages.setup.wavedet.freq):za,2)),2,0,+1);
           %RUTE: can give error in i=1 if za < 0.08*(messages.setup.wavedet.freq)
            paux= modmax(flipud(w(max(1,za-round(QRS_wtol_pantza*(messages.setup.wavedet.freq))):za,2)),2,0,+1);
            if paux,  pantza = za - paux(1)+1;
            else  pantza = []; end
            apantza = w(pantza,2);
            
            %  Protection against waves with cracks
            %paux=modmax(flipud(w(max(1,pantza-round(0.05*(messages.setup.wavedet.freq))):pantza-1,2)),2,umbral2a*dermax,sign(w(pantza,2))); 
            % 0.05 why, when it was changed?
            paux =modmax(flipud(w(max(1,pantza-round(QRS_wtol_pantza_crack*(messages.setup.wavedet.freq))):pantza-1,2)),2,umbral2a*dermax,sign(w(pantza,2)));
            if ~isempty(paux)
                paux = pantza -paux(end);
                pantza = paux;
                apantza = sign(w(paux,2))*max(abs(apantza), abs(w(paux,2)));
            end
            
            if (apantza<0)|(abs(apantza)<umbral1a*dermax)%#ok<OR2>
                za = [];
            elseif isempty(pantza),
                za = [];
            end
            if ~isempty(za), 
                  
               % ind = zerocros(flipud(w(pantza-0.08*(messages.setup.wavedet.freq):pantza-1,1)));
                %%RUTE: can give error in i=1 if pantza < 0.08*(messages.setup.wavedet.freq)
                ind = zerocros(flipud(w(max(1,pantza-round(QRS_wtol_zaa*(messages.setup.wavedet.freq))):pantza-1,1))); %%RUTE:  %07 Abri06
                zaa = pantza-ind;
                if ~isempty(zaa),
                    paux= modmax(flipud(w(max(1,zaa-round(QRS_wtol_pantzaa*(messages.setup.wavedet.freq))):zaa,2)),2,0,-1);
                    if paux,  pantzaa = zaa - paux(1)+1;
                    else  pantzaa = []; end
                    apantzaa = w(pantzaa,2);
                    
                    if (apantzaa>0)|(abs(apantzaa)<umbral2a*dermax), %#ok<OR2>
                        zaa=[];
                    elseif(isempty(pantzaa))
                        zaa=[];
                    end
                end
            end
        end
        
        ind = zerocros(w(pp+1:min(size(w,1),pp+round(QRS_wtol_zp*(messages.setup.wavedet.freq))),2));
        if ~isempty(pp) && ~isempty(ind) %2JUL2011
        zp = pp + ind;
        end
            
        if ~isempty(zp),
            paux= modmax(w(zp:min(size(w,1),zp+round(QRS_wtol_ppostzp*(messages.setup.wavedet.freq))),2),2,0,-1);
            if paux,  ppostzp = zp + paux(1)-1;
            else  ppostzp = []; end
            appostzp = w(ppostzp,2);
            
            %  Protection against waves with cracks
            paux =modmax(w(ppostzp+1:min(size(w,1),ppostzp+round(QRS_wtol_ppostzp_crack*(messages.setup.wavedet.freq))),2),2,umbral2p*dermax,sign(w(ppostzp,2)));  % estamos a ir longe de mais na busca falta protecçao para que nao entre no batimento seguinte!
            if ~isempty(paux)
                paux = ppostzp + paux(end);
                ppostzp = paux;
                appostzp = sign(w(paux,2))*max(abs(appostzp), abs(w(paux,2)));
            end
            
            if (appostzp>0) |(abs(appostzp)<umbral1p*dermax), %#ok<OR2>
                zp=[];
            elseif isempty(ppostzp),
                zp = [];
            end
            if ~isempty(zp),
                ind = zerocros(w(ppostzp+1:min(size(w,1),ppostzp+round(QRS_wtol_zpp*(messages.setup.wavedet.freq))),2));
                zpp = ppostzp + ind;
                if ~isempty(zpp),
                    %ppostzpp = picpost(w(zpp:min(size(w,1),zpp+0.08*(messages.setup.wavedet.freq)),2),zpp);
                    paux= modmax(w(zpp:min(size(w,1),zpp+round(QRS_wtol_ppostzpp*(messages.setup.wavedet.freq))),2),2,0,+1);
                    if paux,  ppostzpp = zpp + paux(1)-1;
                    else  ppostzpp = []; end
                    
                    appostzpp = w(ppostzpp,2);
                    if (appostzpp<0) |(abs(appostzpp)<umbral2p*dermax), %#ok<OR2>
                        zpp=[];
                    elseif(isempty(ppostzpp))
                        zpp=[];
                    end
                end
            end
        end
        
        if (isempty(zaa)) & (qrs-za)>farza*(messages.setup.wavedet.freq), %#ok<AND2> %changed from 0.16 MAY2011
            za = [];   % Si 1º onda m
        end
        if ~isempty(zaa) & (za - zaa)> farzaa*(messages.setup.wavedet.freq), %#ok<AND2> %changed from 0.08 MAY2011
            zaa=[];
        end
        
        if isempty(zpp) & (zp-qrs)>farzp*(messages.setup.wavedet.freq),%#ok<AND2> %changed from 0.13 MAY2011
            zp = [];
        end
        if ~isempty(zpp) & (zpp-zp)>farzpp*(messages.setup.wavedet.freq) %#ok<AND2> %changed from 0.13 MAY2011
            zpp = [];
        end
        
        % wave label asignation
        
        if isempty(za),             % if no waves before qrs
            Q= qrs; R= zp; S=zpp;      % (Q)RS
            if (isempty(zp) && isempty(zpp)),  % no waves before nor after
                Q=[]; S=qrs; end % QS complex
        elseif isempty(zp),         % if some wave before, but not after
            Q=zaa; R=za; S=qrs;        % QR(S) or R(S)
        else                        % if some wave before and some wave after
            if (~isempty(zpp)) && (~isempty(zaa)) % if two waves before and two after
                if (abs(apantzaa)>abs(appostzpp)), % forget the smaller of the extreme waves
                    zpp=[]; 
                else zaa=[];
                end
            end
            % now there are not more than 4 possible waves
            if ~isempty(zpp),                  % if 1 wave before and two after
                if abs(appostzpp)>abs(apantza),   % if the second after greater than the wave before
                    Q = qrs; R = zp;S = zpp;        % (Q)RS
                else                              % if the wave before qrs greater than the secon wave after qrs
                    R = za;                         % R(S)R'
                    S = qrs;
                    Rprima = zp;
                end
            elseif ~isempty(zaa),              % if 2 waves before and one after
                if abs(apantzaa)>abs(appostzp),
                    Q=zaa; R=za; S=qrs;  % QR(S)
                else
                    R=za; S=qrs; Rprima=zp; % R(S)R'
                end
            else                               % if 1 before and one after qrs
                R=za; S=qrs; Rprima=zp; % R(S)R'
            end
            
        end
        
    end 
     
    
    Konset = Kq;   %constants used for onset and offset detection
    Koffset = Ks;
    
    % Identify the first wave in the complex
    
    if ~isempty(Q),
        firstwave = Q;
    elseif ~isempty(R),
        firstwave =R; Konset=Kron;
    else
        firstwave = S;
        Konset = Kron;
    end
    
    % Identify the last wave in the complex
    if ~isempty(Rprima),
        lastwave = Rprima; Koffset = Kroff;
    elseif ~isempty(S),
        lastwave = S;
    else
        lastwave = R;
        Koffset=Kroff;
    end
    
    % The point of departure for detecting the onset and offset of QRS is a peak in the wavelet transform
    % Now we identify which peak must we use for onset and for offset of QRS
    
    if (firstwave == qrs),
        piconset=pa;  
    elseif (firstwave == za),
        piconset=pantza;
    elseif (firstwave == zaa),
        piconset = pantzaa;
    else
        error('Imposible Condition')
    end
    
    
    if (lastwave == qrs),
        picoffset=pp;
    elseif (lastwave == zp),
        picoffset=ppostzp;
    elseif (lastwave == zpp),
        picoffset = ppostzpp;
    else
        error('Imposible Condition')
    end
    
    % Application of derivative criteria: onset of QRS
    
    % Search onset by means of a threshold in the wavelet related to the amplitude of the wavelet at the peak
    %qrson = searchon (piconset, w(piconset-0.08*(messages.setup.wavedet.freq):piconset,2), Konset);
     %%RUTE: can give error in i=1 if piconset < 0.08*(messages.setup.wavedet.freq) 14Mar06
     qrson = searchon (piconset, w(max(1,piconset-round(QRSon_window*(messages.setup.wavedet.freq))):piconset,2), Konset); %14Mar06
    
    % Protection of onsets too far from qrs
    % If qrs onset more than 120 ms before R wave, there is no Q wave.
    if (firstwave==Q)&(Q~=qrs), %#ok<AND2>
        if ((qrs-qrson)>=firstwaverestriction*(messages.setup.wavedet.freq))||(Q-qrson>=Qwaverestriction*(messages.setup.wavedet.freq)),
            firstwave = R; 
            if (firstwave == qrs),
                piconset=pa;  
            elseif (firstwave == za),
                piconset=pantza;
            elseif (firstwave == zaa),
                piconset = pantzaa;    
            end
            % New search
            %qrson = searchon (piconset, w(piconset-0.08*(messages.setup.wavedet.freq):piconset,2), Kron);
            %%RUTE: can give error in i=1 if piconset < 0.08*(messages.setup.wavedet.freq) 14Mar06
            qrson = searchon (piconset, w(max(1,piconset-round(QRSon_window*(messages.setup.wavedet.freq))):piconset,2), Kron);
            Q = [];
        end
    end
    
    % if first wave is R and qrs onset more than firstwaverestriction w before qrs, there is no R wave.
    if (firstwave==R)&(lastwave~=R)& (qrs~=R), %#ok<AND2> %!!!!!!!!!
        if ((qrs-qrson)>firstwaverestriction*(messages.setup.wavedet.freq))||(R-qrson>=Rwaverestriction*(messages.setup.wavedet.freq)), % 0.2 0.1
            
            
            % 19Out09
%             R = []; % 19Out09
%             firstwave = S; % 19Out09
            %check qrs morphology
            if isempty(Rprima)==1 % QS complex (S wave only)
                firstwave = S;
                R=[];
            else % SRprima -> QR complex
                R=Rprima;
                Q=S;
                S=[];
                Rprima=[];
                firstwave = Q;
            end
            if (firstwave == qrs),
                piconset=pa;  
            elseif (firstwave == za),
                piconset=pantza;
            elseif (firstwave == zaa),
                piconset = pantzaa;    
            end
            %qrson = searchon (piconset, w(piconset-0.08*(messages.setup.wavedet.freq):piconset,2), Ks);
            %%RUTE: can give error in i=1 if piconset < 0.08*(messages.setup.wavedet.freq) 14Mar06
            qrson = searchon (piconset, w(max(1,piconset-round(QRSon_window*(messages.setup.wavedet.freq))):piconset,2), Ks); %14Mar06                        
            
        end
    end
    
    
    % Application of derivative criteria: offset of QRS
    
    qrsoff = searchoff (picoffset, w(picoffset:min(size(w,1),picoffset+round(QRSoff_window*(messages.setup.wavedet.freq))),2), Koffset);
    
    % Protection of offsets too far from qrs
    % If Rprima offset more than 120 ms after S wave, there is no Rprima wave.
    if (lastwave==Rprima)&(Rprima~=qrs), %#ok<AND2>
        if ((qrsoff-qrs)>Rprimarestriction1*(messages.setup.wavedet.freq))||(qrsoff-Rprima>Rprimarestriction2*(messages.setup.wavedet.freq)),
            lastwave = S;
            if (lastwave == qrs),
                picoffset=pp;
            elseif (lastwave == zp),
                picoffset=ppostzp;
            elseif (lastwave == zpp),
                picoffset = ppostzpp;
            end
            qrsoff = searchoff (picoffset, w(picoffset:min(size(w,1),picoffset+round(QRSoff_window*(messages.setup.wavedet.freq))),2), Koffset);
            Rprima = [];
        end
    end
    % if S offset more than 200 ms after R wave, there is no S wave.
    if (lastwave==S)&(firstwave~=S)&(S~=qrs), %#ok<AND2>
        if ((qrsoff-qrs)>Srestriction1*(messages.setup.wavedet.freq))||(qrsoff-S>Srestriction2*(messages.setup.wavedet.freq)),
            lastwave = R;
            if (lastwave == qrs),
                picoffset=pp;
            elseif (lastwave == zp),
                picoffset=ppostzp;
            elseif (lastwave == zpp),
                picoffset = ppostzpp;
            end
            qrsoff = searchoff (picoffset, w(picoffset:min(size(w,1),picoffset+round(QRSoff_window*(messages.setup.wavedet.freq))),2), Kroff);
            S = [];
        end
    end
    
    if isempty(piconset), piconset = NaN;  end; % 11.03.09
    if isempty(picoffset), picoffset = NaN; end;% 11.03.09
    
    qrspiconset=[qrspiconset piconset]; %#ok<AGROW>
    qrspicoffset=[qrspicoffset picoffset]; %#ok<AGROW>
    
   
    
    %%%%% OUT09%%%% main positive and negative wave
    if exist('apantzaa','var'), if ~isempty(apantzaa),    WT_extrema(1,i)=apantzaa;  pos_extrema(1,i)=pantzaa;end;end
    if exist('apantza','var'), if ~isempty(apantza),    WT_extrema(2,i)=apantza;  pos_extrema(2,i)=pantza; end;end
    if exist('apa','var'), if ~isempty(apa),    WT_extrema(3,i)=apa;  pos_extrema(3,i)=pa; end;end
    if exist('app','var'), if ~isempty(app),    WT_extrema(4,i)=app; pos_extrema(4,i)=pp; end;end
    if exist('appostzp','var'), if~isempty(appostzp),    WT_extrema(5,i)=appostzp; pos_extrema(5,i)=ppostzp; end;end
    if exist('appostzpp','var'), if ~isempty(appostzpp),    WT_extrema(6,i)=appostzpp; pos_extrema(6,i)=ppostzpp;  end;end
    if exist('zaa','var'), if  ~isempty(zaa),    WT_zc(1,i)=zaa;  end;end
    if exist('za','var'), if ~isempty(za),    WT_zc(2,i)=za;  end;end
    if exist('qrs','var'), if ~isempty(qrs),    WT_zc(3,i)=qrs; end;end
    if exist('zp','var'), if ~isempty(zp),    WT_zc(4,i)=zp;  end;end
    if exist('zpp','var'), if ~isempty(zpp),    WT_zc(5,i)=zpp;  end;end

    
    
    
    %%%%%%%% RUTE01MAR2012
    aux=[];
    if ~isempty(R)
        aux=find(WT_zc(:,i)==R,1);
    end
    if ~isempty(Q)
        aux=[aux find(WT_zc(:,i)==Q,1)]; %#ok<AGROW>
    end
    if ~isempty(S)
        aux=[aux find(WT_zc(:,i)==S,1)]; %#ok<AGROW>
    end
    if ~isempty(Rprima)
        aux=[aux find(WT_zc(:,i)==Rprima,1) aux]; %#ok<AGROW>
    end
    
    auxaux=NaN(5,1);
    auxaux(aux)=WT_zc(aux,i);
    WT_zc(:,i)=auxaux;
    %%%%%%%%%%%%%%%
    
    A=find(~isnan(WT_zc(:,i)));
    auxaux= WT_extrema(A,i)-WT_extrema(A+1,i);
    [M,iM]=max(auxaux); %#ok<ASGLU>
    [m,im]=min(auxaux); %#ok<ASGLU>
    pos.QRSmainpos(i)= WT_zc(A(iM),i);
    pos.QRSmaininv(i)= WT_zc(A(im),i);
    
    % maximum slope locations for main wave % 2AGO2011  
    if ~isempty(pa)
        pos.QRSpa(i)=pa;
    else
        pos.QRSpa(i) = NaN;
    end
    if ~isempty(pp)        
        pos.QRSpp(i)=pp;
    else
        pos.QRSpp(i) = NaN;
    end
    %%%%%%%%%%
    
    %  Filling the structs of positions
    if isempty(qrson), qrson = NaN;  end;
    if isempty(qrsoff), qrsoff = NaN; end;
    if isempty(R), R=NaN; end;
    if isempty(S), S=NaN; end;
    if isempty(Q), Q=NaN; end;
    if isempty(Rprima), Rprima=NaN; end;
    
    pos.QRSon(i)= qrson;
    pos.Q(i) = Q;
    pos.R(i) = R;
    pos.S(i) = S;
    pos.Rprima(i) = Rprima;
    pos.QRSoff(i) = qrsoff;
    pos.qrs(i) = qrs;  
        
end

if ~isempty(intervalo)
    position.QRSon(intervalo(1):intervalo(2)) = pos.QRSon+samp(1)-1;
    position.Q(intervalo(1):intervalo(2)) = pos.Q+samp(1)-1;
    position.R(intervalo(1):intervalo(2)) = pos.R+samp(1)-1;
    position.S(intervalo(1):intervalo(2)) = pos.S+samp(1)-1;
    position.Rprima(intervalo(1):intervalo(2)) = pos.Rprima+samp(1)-1;
    position.QRSoff(intervalo(1):intervalo(2)) = pos.QRSoff+samp(1)-1;
    position.QRSmainpos(intervalo(1):intervalo(2))= pos.QRSmainpos+samp(1)-1; %28OUT09
    position.QRSmaininv(intervalo(1):intervalo(2))= pos.QRSmaininv+samp(1)-1; %28OUT09
    position.QRSpa(intervalo(1):intervalo(2))=pos.QRSpa+samp(1)-1;% 2AGO2011
    position.QRSpp(intervalo(1):intervalo(2))=pos.QRSpp+samp(1)-1;% 2AGO2011
end