function position = wavedet(sigdir,headir,matdir,ecgnr,ft,anot,lead,t,qrs_flag,anot_fmt,aname,dirann)
%function position = mainwav(sigdir, headir, matdir, ecgnr, ft, anot, lead, t)
%Significant ECG points detector based on wavelet approach %Syntaxis:
%Input Parameters:

%	sigdir: directory of the original signal
%	headir: directory of the header file
%	matdir: output directory
%	ecgnr: string with the name of the record to analyze
% 	ft: format file: 0 for MIT header 1 for LUND header 2 for a mat file
% 	without header.
%	anot: name of the annotation output file
%	lead: scalar with the lead number to analyze (1,2, ...) (only one lead)
%	t=[tbegin tend]: time vector with initial and sample to analyze
%    qrs_flag: External (1) or internal (0) QRS detector
%    anot_fmt: in case of external QRS detector, format of annotation file:
%    aname: annotation file name with QRS marks
%    dirann: directory with external annotation file
%Output Parameters:
%	position: struct vector with the detected points


% The analysis loop has an overlapping structure.  The number of samples
% processed is large and there are an overlap at the beginning and another
% at the end.  The first one is necessary because: a) the first l4-1 samples
% are always incorrectly filtered b) To align the signals, as the filters
% have different lengths, we must discard some samples. c) The algorithms
% must have some possibility of turning back.
% The second overlap is mainly because of the alignment, as the filters have
% different delays.
% Juan Pablo Martínez Cortés
% 1/3/2006 option ft=2 was added from wavedetplus.
%%%%%alteraçoes para estudar a importancia do limiar de regularidade %%%%%%%
%%%%%% Rute 04/09/02
global regularity
regularity.time=[];
regularity.amp2=[];
regularity.amp3=[];
regularity.amp4=[];
regularity.alphalinha=[];
regularity.alpha1=[];
regularity.alpha2=[];
regularity.alpha3=[];
regularity.alpha4=[];  %%%%% para incluir a escala 5 %% Rute 22/05/02
regularity.thalpha1=[];
regularity.thalpha2=[];
regularity.thalpha3=[];
regularity.th2=[];
regularity.th3=[];
regularity.th4=[];

regularity.thlinha=[];
regularity.escolhidos=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off
ultimo_anot=1;
timeoffset=0;
leadsfile = [];
nsamp = 2^16;   % Number of samples per excerpt for the WT

if ft==0   % MIT format files
    if isunix,          %%%%%%% Rute 13/03/02
        sep = '/';
    else sep = '\';     %%%%%%% Rute 13/03/02
    end                 %%%%%%% Rute 13/03/02

    heasig = readheader([headir sep ecgnr '.hea']);

    %%%%%%%%%%%%%%  ANA  PARA BIOPAC %%%%%%%%%%%%%%%%%%%%%%
    if isfield(heasig,'spf')
        heasig.freq=heasig.freq*heasig.spf(lead);
        heasig.nsamp=heasig.nsamp*heasig.spf(lead);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (heasig.fmt(lead)==16) || (heasig.fmt(lead)==212) || (heasig.fmt(lead)==61)
        formato = num2str(heasig.fmt(lead));
    else
        error('This format is not supported by the program')
    end

    for kk=1:heasig.nsig
        leadsfile(kk)=strcmp(heasig.fname(lead,:),heasig.fname(kk,:));
    end
    leadsfile = cumsum(leadsfile);
    leadinfile = leadsfile(lead);  % this is necessary when there are leads in different files.   e.g. ptbdb database  JP 2006
    nsiginfile = leadsfile(end);

    if strcmp(formato,'16')||strcmp(formato,'61'),
        % fid = fopen([sigdir ecgnr '.dat'],'rb');
        % if fid == -1,
        fid = fopen([sigdir heasig.fname(lead,:)],'rb');
        % end
        fseek(fid,0,-1);  % Rewind the file
        if strcmp(heasig.fname(lead,1),'_'),  % Siemens card recordings with MIT-type header
            timeoffset=512;
            fseek (fid, timeoffset,-1); % offset
        end
    end
elseif ft==1  % LUND format files
    hdsig=gethdsig(sigdir,ecgnr); % Reading header information
    heasig = hdsig2heasig(hdsig);       
    freq = heasig.freq;
elseif ft==2  % data in a Matlab file
    aux=[sigdir ecgnr];
    load ([aux '.mat'])
    heasig.nsamp=length(sinal);
    heasig.freq=fa;
    freq = heasig.freq;
    if exist('gain','var')
        heasig.gain=gain;
    else
        heasig.gain=200;  % When heasig.gain = 0 => default 200
    end
else
    error('Format ft should be 0, 1 or 2')
end
heasig.gain(heasig.gain==0)=200;  % When heasig.gain = 0 => default 200

if nargin >=8,
    t=[max(t(1),1) min(heasig.nsamp,t(2))];
else
    qrs_flag=0;% default value was missing 01/04/02 Rute
    t=[1 heasig.nsamp];
end

% Reading of external QRS annotation file if qrs_flag
if qrs_flag==1     % The QRS fiducial point is read from external annotator
    switch (anot_fmt);
        case 0        % MIT annotation file
            if exist([dirann ecgnr '.' aname],'file')
                s=readannot([dirann ecgnr '.' aname],t);
                s=isqrs(s,heasig,t);
                ext_anot=s.time';
            else error('QRS annotation file not found');
            end
        case 1        % mat file
            if exist([dirann aname '.mat'],'file')
                s=load ([dirann aname]);
                kk=getfield(s);
                eval(['ext_anot=s.' kk ';']);
            else error('QRS annotation file not found');
            end
        otherwise
            error('Bad annotation format for external QRS annotation file');
    end
end

inisamp = t(1);
endsamp = 0;
timeqrs = [];
lastqrs = 1;


if qrs_flag==0 % estimated maximum number of beats
    maxlength = (t(2)-t(1))/heasig.freq*2;
else maxlength=length(ext_anot);
end

%if qrs_flag==1,
%  sel=find(ext_anot<heasig.freq);
%  if ~isempty(sel)
%      ext_anot(sel)=[];
%      maxlength= length(ext_anot);
%  end
%end

nanvec = nan*ones(1,maxlength);
position=struct('Pon',nanvec,'P',nanvec,'Poff',nanvec,'QRSon',nanvec,...
    'Q',nanvec,'R',nanvec,'Fiducial',nanvec,'qrs',zeros(1,maxlength),...
    'Rprima',nanvec,'S',nanvec,'QRSoff',nanvec,'Ton',nanvec,'T',nanvec,...
    'Tprima', nanvec,'Toff',nanvec,'Ttipo',nanvec);

if qrs_flag==1, position.qrs=ext_anot;  end


numlatdet = 0;      % Number of detected beats until now

[q1,q2,q3,q4,q5] = qspfilt5(heasig.freq);  % WT equivalent filters %%%%% para incluir a escala 5 %% Rute 22/05/02
l1=length(q1);l2=length(q2);l3=length(q3);l4=length(q4);
d1=floor((l1-1)/2);d2=floor((l2-1)/2);
d3=floor((l3-1)/2);d4=floor((l4-1)/2);
l5=length(q5);d5=floor((l5-1)/2); %%%%% para incluir a escala 5 %% Rute 22/05/02
begoverlap = l5+2*heasig.freq;    % l5 + 2 sec. %%%%% passa a ser em relaçao a escala 5 %% Rute 22/05/02
endoverlap = d5;                                %%%%% passa a ser em relaçao a escala 5 %% Rute 22/05/02
% The two seconds overlap makes sure that the last/first beat
% are completely within the excerpt processed

samp=1; %freq;    %Initialization of samp (1 sec)

while ((endsamp+1) < t(2))
    firstnewsamp = samp(end);   % last sample analyzed
    endsamp = min(inisamp + nsamp -1,t(2));

    if ft==0 % MIT format
        if strcmp(formato,'212'),  % !!!!! Different formats (heasig)
            sig = rdsign212([sigdir ecgnr '.dat'],nsiginfile,inisamp, endsamp);  % I changed heasig.nsig by nsiginfile  ---> see above JP2006
        elseif strcmp(formato,'16')|strcmp(formato,'61'),

            %%%%%%%%%%%%%%  ANA  PARA BIOPAC %%%%%%%%%%%%%%%%%%%%%%
            if isfield(heasig,'spf')
                if(~isempty(heasig.spf) & find(heasig.spf ~= 1) )
                    agrupacion = sum(heasig.spf);
                    %         ecg = zeros(heasig.nsig,heasig.nsamp*(max(heasig.spf)));
                end
                fid = fopen([sigdir ecgnr '.ecg'],'rb');
                %      fid = fopen([prmBrowse.ecgdir prmBrowse.ecgnr],'rb');
                if(length(find(heasig.fmt == 16)) == heasig.nsig)
                    %fseek(fid,agrupacion*(inisamp-1)*2,'bof');
                    fseek(fid,agrupacion*round((inisamp-1)/heasig.spf(lead))*2,'bof');
                    x = fread(fid,[agrupacion (endsamp-inisamp+1)/heasig.spf(lead)], 'int16');
                    ini_pos_lead = sum(heasig.spf(1:lead))-heasig.spf(lead)+1;
                    end_pos_lead = sum(heasig.spf(1:lead));
                    sig = reshape(x(ini_pos_lead:end_pos_lead,:),[],1);
                    % Para no tener que cambiar el codigo
                    sig=[zeros(length(sig),leadinfile-1) sig];
                    clear('x')
                end
                fclose('all');
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            else
                fseek(fid, timeoffset+2*(inisamp-1)*nsiginfile, -1); % Locate the pointer
                sig = fread(fid,[nsiginfile endsamp-inisamp+1] ,'int16')';
            end
            if (strcmp(computer,'SOL2')&strcmp(formato,'16'))|(~strcmp(computer,'SOL2')&strcmp(formato,'61')),
                sig=swap16(sig);
            end
        end
        sig = sig(:,leadinfile);   % Only selected lead.  --> JP 2006  I changed lead by leainfile... see at the top

        if exist('heasig','var') && isfield(heasig,'units')
            if upper(heasig.units(lead))=='MV'
                sig = (sig - heasig.adczero(lead))/heasig.gain(lead)*1e3; % conversion to microV
            elseif upper(heasig.units(lead))=='V'
                sig = (sig - heasig.adczero(lead))/heasig.gain(lead)*1e6; % conversion to microV
            else
                sig = (sig - heasig.adczero(lead))/heasig.gain(lead);
            end
        else
            sig = (sig - heasig.adczero(lead))/heasig.gain(lead)*1e3; % conversion to microV

        end

        %         sig = (sig - heasig.adczero(lead))/heasig.gain(lead)*1e3;  % conversion to microV
    elseif ft==1  % LUND format
        sig=getsig(sigdir,ecgnr,[inisamp,endsamp],lead)';
        %sig=sig/5;   % Mejor en uV JP
    elseif ft==2  % data in a Matlab file
        sig=sinal(lead,inisamp:endsamp)';
    end

    wt = wavt5(sig',q1,q2,q3,q4,q5);  % WT equivalent filtering % 22/05/02 Rute including scale 5
    wt=wt(l5:end,1:5);            % Remove "incorrect" samples
    samp = (inisamp + l5 -1):endsamp-d5;
    % First l5-1 samples are not correctly filtered (border effect)
    % Last d5 samples are discarded in order to alineate all the
    % filtered signals, taking into acount the filter delays

    % Synchronizing filtered signals at different scales
    w = zeros(size(wt,1)-d5,5);
    w(:,5) = wt(d5+1:end,5);
    w(:,4) = wt(d4+1:end+d4-d5,4);
    w(:,3) = wt(d3+1:end+d3-d5,3);
    w(:,2) = wt(d2+1:end+d2-d5,2);
    w(:,1) = wt(d1+1:end+d1-d5,1);
    sig=sig(l5:end-d5); % piece of signal processed in this iteration

    clear wt;
    % Initialization of the thresholds for QRS detection (eps)
    eps= 0.5*sqrt(mean(w.^2));
    eps(4) = eps(4)*2;

    if qrs_flag==0  % The QRS fiducial point is calculated with wavedet
        %fiducial_altera_novo;  % uso de Criterios alternativos para selecao de candidatos a QRS % Rute 4/9/02
        fiducial; % Detection of fiducial point % 22/05/02 Rute including scale 5 % 07/06/02 Rute exclui escala5
    else           % The QRS fiducial point is read from external annotator
        first = max(firstnewsamp-ceil(heasig.freq), samp(1)-1+ceil(heasig.freq*0.050));
        % first is calculated according with the first lines of fiducial JP

        sel=find(position.qrs>=first & position.qrs<=samp(end));
        time=position.qrs(sel)-samp(1)+1;


        % If the P-wave is not in the present segment, discard the beat, because the
        % the beat should have been detected in the last segment.
        if abs(samp(time(1)) - lastqrs)< 0.275*heasig.freq,
            time(1)=[];  sel(1)=[];              % In order not to repeat beats
        end

        if (length(sig)-time(end) < 0.45*heasig.freq),
            time(end)=[]; sel(end)=[];
        end
        % If last beat's T-wave is not in 'sig', the beat is detected in next one

        if time(1)-0.35*heasig.freq<=0, time(1)=[]; sel(1)=[]; end

        timeqrs = [timeqrs samp(time)];  % QRS times

        lastqrs = samp(time(end));

        intervalo = [sel(1) sel(end)];       %%%%%%
        % numeration of the beats processed in this segment
        %keyboard
    end
    if ~isempty(time)
        %% Individual wave detection  %%
        qrswave;
        %%%%%%%%%%%%%%%% caso em que rr tem dimensao inferior ao necessario
        if length(time)>1  %%%%%%%%%%%%%%%%introduzido a 17/05/02 Rute
            twave,  %%%%%%%% 07/06/0 Rute including scale 5
        end %%%%%%%%%%%%%%%%introduzido a 17/05/02 Rute
        pwave;
        % Last annotated position
        if ~isnan(position.Toff(intervalo(2))),
            ultimo_anot = position.Toff(intervalo(2));
        elseif ~isnan(position.Tprima(intervalo(2))),
            ultimo_anot = position.Tprima(intervalo(2));
        elseif ~isnan(position.T(intervalo(2))),
            ultimo_anot = position.T(intervalo(2));
        elseif ~isnan(position.QRSoff(intervalo(2))),
            ultimo_anot = position.QRSoff(intervalo(2));
        else
            ultimo_anot = position.qrs(intervalo(2));
        end
    end

    inisamp = endsamp+2-endoverlap-begoverlap;
end


% Remove void annotations at the end
aux = find (position.qrs == 0);
position.Pon(aux)=[];
position.P(aux)=[];
position.Poff(aux)=[];
position.QRSon(aux)=[];
position.Q(aux)=[];
position.R(aux)=[];
position.Fiducial(aux)=[];
position.qrs(aux)=[];
position.Rprima(aux)=[];
position.S(aux)=[];
position.QRSoff(aux)=[];
position.Ton(aux)=[];
position.T(aux)=[];
position.Tprima(aux)=[];
position.Toff(aux)=[];
position.Ttipo(aux)=[];


annstruct = pos2ann(position);
if ~isempty(position.qrs)
    writeannot ([matdir ecgnr '.' anot 's' num2str(lead)], annstruct);   % mexfile, faster!
end

%escribeanot(fann,position,anot,lead, ecgnr, freq, matdir, 'anot');

fclose('all');
clear global regularity