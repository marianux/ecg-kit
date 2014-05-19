function     [position,timeqrs,lastqrs,intervalo,numlatdet,time,messages]= fiducialf(position,firstnewsamp,samp,heasig,w,eps,nlead,lastqrs,timeqrs,numlatdet,ultimo_anot,messages)
% [position,timeqrs,lastqrs,intervalo,numlatdet,time,messages]= fiducialf(position,firstnewsamp,samp,heasig,w,eps,nlead,lastqrs,timeqrs,numlatdet,messages)
%
% This script performs the QRS fiducial point detection of ECG
% In this version the thresholds eps(x) are not beat to beat actualized.
%
%Input Parameters:
%   position: struct vector with the detected points
%   firstnewsamp: last sample analyzed in the previous excerpt
%   samp: samples included in the current excerpt (borders excluded)
%   heasig: struct vector with header information
%   w: matrix with WT scales 1 to 5
%   eps:  thresholds
%   nlead: lead to be delineated
%   lastqrs: last detecter QRS time (shoud be timeqrs(end))
%   timeqrs:  QRS times
%   numlatdet: Number of detected beats until now
%   ultimo_anot: last annotation in the previous segment

%Output Parameters:
%   intervalo: numeration of the beats processed in this segment
%   time: beats processed in this segment
%   Actualized parameters: position,timeqrs,lastqrs,numlatdet,timeqrs
%
% Last update: Rute Almeida 07FEB2012
%
% Designed for MATLAB Version R12; tested with MATLAB Version R13
if nargin<12 || ~isfield(messages.setup,'wavedet')
    messages.setup.wavedet=[];
end
if isfield(messages.setup.wavedet,'auxff') && messages.setup.wavedet.auxff==1
    ff=1;
    if isfield(messages.setup.wavedet,'aux_aa')
        aa=messages.setup.wavedet.aux_aa;
    else
        aa=[1 1500];
    end
else
    ff=0;
end
if isempty(heasig)
    
end
if ff==1
    figure
    plot(w(:,1),'k')
    fig1=gca;
    set(gca,'xlim',aa)
    hold on
    ylabel('W_{2^1}')
    figure
    plot(messages.sig,'k')
    fig2=gca;
    hold on
    set(gca,'xlim',aa)
    ylabel('ECG')
    figure
    plot(w(:,3),'k')
    fig3=gca;
    set(gca,'xlim',aa)
    hold on
    ylabel('W_{2^3}')    
    figure
    subplot(4,1,1)
    plot(w(:,4),'k')
    ylabel('W_{2^4}')
    hold on
    figWT4=gca;
    set(gca,'xlim',aa)
    subplot(4,1,2)
    plot(w(:,3),'k')
    ylabel('W_{2^3}')
    hold on
    figWT3=gca;
    set(gca,'xlim',aa)
    subplot(4,1,3)
    plot(w(:,2),'k')
    hold on
    ylabel('W_{2^2}')
    figWT2=gca;
    set(gca,'xlim',aa)
    subplot(4,1,4)
    plot(w(:,1),'k')
    hold on
    ylabel('W_{2^1}')
    figWT1=gca;
    set(gca,'xlim',aa)
    
    figure
    subplot(4,1,1)
    plot(messages.sig,'k')
    ylabel('ECG')
    hold on
%     figECG=gca;
    set(gca,'xlim',aa)
    subplot(4,1,2)
    plot(w(:,3),'k')
    ylabel('W_{2^3}')
    hold on
    figECGWT3=gca;
    set(gca,'xlim',aa)
    subplot(4,1,3)
    plot(w(:,2),'k')
    hold on
    ylabel('W_{2^2}')
    figECGWT2=gca;
    set(gca,'xlim',aa)
    subplot(4,1,4)
    plot(w(:,1),'k')
    hold on
    ylabel('W_{2^1}')
    figECGWT1=gca;
    set(gca,'xlim',aa)
    
end
if~isfield(messages.setup.wavedet,'firstmin')
    messages.setup.wavedet.firstmin=0.050; % minimum time in sec to discart at the beggining
end
if ~isfield(messages.setup.wavedet,'nghbhd')
    messages.setup.wavedet.nghbhd=0.025; % neighbourwood for maximum across scales in sec
end
if ~isfield(messages.setup.wavedet,'peakcriteria')
    messages.setup.wavedet.peakcriteria=1.2; % If a modulus is more than peakcriteria times greater than others in the sane scale the greatest is chosen
end
if ~isfield(messages.setup.wavedet,'thmax')
    messages.setup.wavedet.thmax=5+log(5); %maximum value for log(abs(W_3)) in a maximun line
end
if ~isfield(messages.setup.wavedet,'th_schb')
    messages.setup.wavedet.th_schb=5+log(5); %maximum value for log(abs(W_3)) in a maximun line
end
if ~isfield(messages.setup.wavedet,'thfraction')
    messages.setup.wavedet.thfraction=0.5;   % tolerance below the average log(abs(W_3)) in a maximun line
end
if ~isfield(messages.setup.wavedet,'intvlthr2')
    messages.setup.wavedet.intvlthr2=0.15;%threshold interval to consider isolated maximum lines
end
if ~isfield(messages.setup.wavedet,'intvlthr1')
    messages.setup.wavedet.intvlthr1=0.12;%threshold interval for considering redundancy
end
if ~isfield(messages.setup.wavedet,'intvlthr1_2')
    messages.setup.wavedet.intvlthr1_2=1.5;%reduction on threshold interval for considering redundancy if it subsists
end
if ~isfield(messages.setup.wavedet,'pictime')
    messages.setup.wavedet.pictime=0.1;%length of the search window for associated maximumlines in extra verification step
end
if ~isfield(messages.setup.wavedet,'timelapthr')
    messages.setup.wavedet.timelapthr=0.125;% maximum interval between two assiciated maximumlines in sec
end
if ~isfield(messages.setup.wavedet,'refrper')
    messages.setup.wavedet.refrper=0.275;% Refractary period after a QRS detection in sec
end
if ~isfield(messages.setup.wavedet,'refrper_reduc')% 22NOV2011
    messages.setup.wavedet.refrper_reduc=0;% 22NOV2011
end
if ~isfield(messages.setup.wavedet,'rrmaxmaxlim')
    messages.setup.wavedet.rrmaxmaxlim=0.3; % search back criteria: searchback if RR>wavedet.rrmaxmaxlim sec
end
if ~isfield(messages.setup.wavedet,'rrmaxlim')
    messages.setup.wavedet.rrmaxlim=1.5;% search back criteria: searchback if RR>wavedet.rrmaxlim * median(RR)
end
if ~isfield(messages.setup.wavedet,'rrmaxlimsec')
    messages.setup.wavedet.rrmaxlimsec=1.5;% search back criteria: searchback if RR>wavedet.rrmaxlimsec
end
if ~isfield(messages.setup.wavedet,'rrmax_rbt')
    messages.setup.wavedet.rrmax_rbt=0.1; % search back criteria only if >wavedet.rrmax_rbt sec o
end
if ~isfield(messages.setup.wavedet,'gapthrs')
    messages.setup.wavedet.gapthrs=2; % fraction of mean(rr) to decide if there are a final "gap"
end
if ~isfield(messages.setup.wavedet,'pwaveout')
    messages.setup.wavedet.pwaveout=0.2;% mimimum time in sec before each beat in order to ensure a complete beat
end
if ~isfield(messages.setup.wavedet,'twaveout')
    messages.setup.wavedet.twaveout=0.45; % mimimum time in sec after each beat in order to ensure a complete beat
end
if ~isfield(messages.setup.wavedet,'refrperdef')
    messages.setup.wavedet.refrperdef=messages.setup.wavedet.refrper;%
end
if messages.setup.wavedet.refrper_reduc>=1 % 22NOV2011
    messages.setup.wavedet.refrper_reduc=0.9;% 22NOV2011
end

refrper=messages.setup.wavedet.refrper(end);
peakcriteria=messages.setup.wavedet.peakcriteria;
%global regularity
intervalo=[];
first = firstnewsamp-round(samp(1))+1;
first = max(ceil(messages.setup.wavedet.freq*messages.setup.wavedet.firstmin), first - ceil(1*messages.setup.wavedet.freq));

n = modmax(w(:,4),first,eps(4,nlead),0); % Maximum moduli at scale 4 % Rute 02.Dec.04

if ~isempty(n),
    signo = sign(w(n,4))';             % Maximum or minimum?
    neighb = ceil(messages.setup.wavedet.freq*messages.setup.wavedet.nghbhd);  % Definition of neighborhood 25 ms.
    %it makes sence to verify if all m(4,:) diffre at least neighb....
    m = zeros(4,length(n));
    m(4,:) = n';                       % Positions of maximum moduli
    for k = 1:size(m,2),               % For each maximum at scale 4
        n3 = modmax(w(max(1,n(k)-neighb):min(n(k)+neighb,size(w,1)),3),2,eps(3,nlead),signo(k));% Rute 02.Dec.04
        n3 =n(k)-neighb-1+n3;            % Check the neighborhood at scale 3
        num = length(n3);
        if num>0,
            if num==1,                     % If only one neighbour maximum
                m(3,k)= n3;                   % at scale 3
            elseif num>1                   % If more than one
                if ~isempty(find(max(abs(w(n3,3))./(abs(w(n3,3)))<peakcriteria))==1),
                    % If a modulus is more than peakcriteria times greater than others
                    [aux1,ind]=max(abs(w(n3,3))); %#ok<ASGLU> % the greatest is chosen
                    m(3,k)= n3(ind);
                else                          % Choose the one with minimum distance
                    [aux1,ind]=min(abs(m(4,k)-n3)); %#ok<ASGLU>
                    m(3,k)= n3(ind);
                end
            end
        end
    end
    ind = find(m(3,:)==0);            % Discard all maximum lines with no
     if ff==1
       fig4WT = axes(figWT4);
       plot(fig4WT,m(4,:),w(m(4,:),4),'ob')
       plot(fig4WT,[aa(1) aa(2)],[eps(4,nlead) eps(4,nlead)],':b')
       plot(fig4WT,[aa(1) aa(2)],-[eps(4,nlead) eps(4,nlead)],':b')
       if ~isempty(ind)
           plot(fig4WT,m(4,ind),w(m(4,ind),4),'xr')
           plot(fig4WT,[m(4,ind(1))-neighb m(4,ind(1))-neighb],[-5000 5000],':r')
           plot(fig4WT,[m(4,ind(1))+neighb m(4,ind(1))+neighb],[-5000 5000],':r')
           fig3WT = axes(figWT3);
           plot(fig3WT,[m(4,ind(1))-neighb m(4,ind(1))-neighb],[-5000 5000],':r')
           plot(fig3WT,[m(4,ind(1))+neighb m(3,ind(1))+neighb],[-5000 5000],':r')
       end
       if ~exist(fig3WT,'var')
           fig3WT = axes(figWT3);
       end
       plot(fig3WT,m(3,m(3,:)~=0),w(m(3,m(3,:)~=0),3),'ob')
       plot(fig3WT,[aa(1) aa(2)],[eps(3,nlead) eps(3,nlead)],':b')
       plot(fig3WT,[aa(1) aa(2)],-[eps(3,nlead) eps(3,nlead)],':b')
     end
    m(:,ind) = [];                    % associated maximum at scale 3
    signo(ind)= [];

    % Search for maximum modula in the neighborhood at scale 2
    for k = 1:size(m,2),
        n2 =modmax(w(max(1,m(3,k)-neighb):min(m(3,k)+neighb,size(w,1)),2),2,eps(2,nlead),signo(k));% Rute 02.Dec.04
        n2 =m(3,k)-neighb-1+n2;
        num = length(n2);
        if num>0,
            if num==1, 			     %  If only one neighbour maximum
                m(2,k)= n2;  		     %  at scale 2
            elseif num>1                     % If more than one
                if ~isempty(find(max(abs(w(n2,2))./(abs(w(n2,2)))<peakcriteria))==1),
                    [aux1,ind]=max(abs(w(n2,2))); %#ok<ASGLU> % greatest modulus
                    m(2,k)= n2(ind);
                else                            % shortest distance
                    [aux1,ind]=min(abs(m(3,k)-n2)); %#ok<ASGLU>
                    m(2,k)= n2(ind);
                end
            end
        end
    end
    
    ind = find(m(2,:)==0);		     % Discard all maximum lines with no
   
    if ff==1
%         axes(figWT3);
        if ~isempty(ind)
            plot(fig3WT,[m(3,ind(1))-neighb m(3,ind(1))-neighb],[-5000 5000],':r')
            plot(fig3WT,[m(3,ind(1))+neighb m(3,ind(1))+neighb],[-5000 5000],':r')
            plot(fig3WT,m(3,ind),w(m(3,ind),3),'xr')
            fig2WT = axes(figWT2);
            plot(fig2WT,[m(3,ind(1))-neighb m(3,ind(1))-neighb],[-5000 5000],':r')
            plot(fig2WT,[m(3,ind(1))+neighb m(3,ind(1))+neighb],[-5000 5000],':r')
            
        end
        if ~exist(fig2WT,'var')
            fig2WT = axes(figWT2);
        end
        plot(fig2WT,m(2,m(2,:)~=0),w(m(2,m(2,:)~=0),2),'ob')
        plot(fig2WT,[aa(1) aa(2)],[eps(2,nlead) eps(2,nlead)],':b')
       plot(fig2WT,[aa(1) aa(2)],-[eps(2,nlead) eps(2,nlead)],':b')
    end
    
    m(:,ind) = [];                      % associated maximum at scale 2
    signo(ind)=[];

    for k = 1:size(m,2),
        n1 = modmax(w(max(1,m(2,k)-neighb):min(m(2,k)+neighb,size(w,1)),1),2,eps(1,nlead),signo(k)); % Rute 02.Dec.04
        n1 =m(2,k)-neighb-1+n1;
        n1 = n1(n1 > 0 & n1 <= size(w,1) );
        num = length(n1);
        if num>0,
            if num==1, 		      %  If only one neighbour maximum
                m(1,k)= n1;                     %  at scale 1
            elseif num>1                     % If more than one
                if ~isempty(find(max(abs(w(n1,1))./(abs(w(n1,1)))<peakcriteria))==1),
                    [aux1,ind]=max(abs(w(n1,1))); %#ok<ASGLU> % greatest modulus
                    m(1,k)= n1(ind);
                else                             % shortest distance
                    [aux1,ind]=min(abs(m(2,k)-n1)); %#ok<ASGLU>
                    m(1,k)= n1(ind);
                end
            end
        end
    end
    ind = find(m(1,:)==0);		     %Discard all maximum lines with no
  
       
    if ff==1
%         axes(figWT2);
        if ~isempty(ind)
            plot(fig2WT,m(2,ind),w(m(2,ind),2),'xr')
        end
        fig1WT = axes(figWT1);
        plot(fig1WT,m(1,m(1,:)~=0),w(m(1,m(1,:)~=0),1),'ob')
        plot(fig1WT,[aa(1) aa(2)],[eps(1,nlead) eps(1,nlead)],':b')
       plot(fig1WT,[aa(1) aa(2)],-[eps(1,nlead) eps(1,nlead)],':b')
    end
    
    m(:,ind) = [];  		     % associated maximum at scale 1
    signo(ind)=[];
    %%%  Regularity Exponent Validation %%%%%%%%%
    % alpha is proportional to log(a3(nk3))-log(a1(nk1))
    alpha = log(abs(w(m(3,:),3))); % !!!!
    %     alphalinha = log(abs(w(m(3,:),3))) - log (abs(w(m(1,:),1))); % !!! (1) % cometario na versao original%%%%%% Rute 24/04/02
    %
    %     %%%% SO
    %     alpha1 = log(abs(w(m(2,:),2))./abs(w(m(1,:),1)));
    %     alpha2 = log(abs(w(m(3,:),3))./abs(w(m(2,:),2)));
    %     alpha3 = log(abs(w(m(4,:),4))./abs(w(m(3,:),3)));
    %     %%%% SO: you are right alpha=alpha1+alpha2 is calculated acording to (1)
    %     %%%% SO: por ahorro computacional se puede eliminar el log, y en el
    %     %%%% SO: valor del threshold si queremos.
    %
    th=min(messages.setup.wavedet.thmax,mean(alpha)-messages.setup.wavedet.thfraction);
    %    thlinha=min(5+log(5),mean(alphalinha)-.5); %%%%%% Rute 24/04/02
    ind = find(alpha<=th); % always empty?????!????1
    % ind=find(alpha3>=1 || alpha1<=0.9 || alpha2<=0.9);% elimination of large T waves
    if ff==1
       figax1 = axes(fig1);
        plot(figax1,m(1,:),w(m(1,:)),'ob')
        legend(figax1,{'WT_1','mm lines'},'Location','NorthEastOutside')
        plot(m(1,ind),w(m(1,ind)),'xr')
        if ~isempty(ind)
            ll={'mm lines','rejected mm lines'};
        else
            ll={'mm lines'};
        end
        figWT1ECG = axes(figECGWT1);
        plot(figWT1ECG,m(1,:),w(m(1,:),1),'ob')
        plot(figWT1ECG,m(1,ind),w(m(1,ind),1),'xr')
        figWT2ECG = axes(figECGWT2);
        plot(figWT2ECG,m(2,:),w(m(2,:),2),'ob')
        plot(figWT2ECG,m(2,ind),w(m(2,ind),2),'xr')
        figWT3ECG = axes(figECGWT3);
        plot(figWT3ECG,m(3,:),w(m(3,:),3),'ob')
        plot(figWT3ECG,m(3,ind),w(m(3,ind),3),'xr')
        
        figax3 = axes(fig3);
        plot(figax3,m(3,:),w(m(3,:),3),'ob')
        legend(figax3,{'WT_3','mm lines'},'Location','NorthEastOutside')
        plot(m(3,ind),w(m(3,ind),3),'xr')
%         axes(fig1);
        legend(figax1,['WT_1',ll],'Location','NorthEastOutside')
%         axes(fig3);
        legend(figax3,['WT_3',ll],'Location','NorthEastOutside')
        
        figure;
        hs = subplot(3,1,1);
        plot(hs,w(:,3),'k')
        ylabel(hs,'W_{2^3}')
        hold on
        plot(hs,m(3,:),w(m(3,:),3),'ob')
        plot(hs,m(3,ind),w(m(3,ind),3),'xr')
        set(hs,'xlim',aa)
        hs3 = subplot(3,1,3);
        plot(hs3,messages.sig,'k')
        hold on
        set(hs3,'xlim',aa)
        ylabel(hs3,'ECG')
        set(hs3,'xlim',aa)
        hs2 = subplot(3,1,2);
        plot(hs2,m(3,:),alpha,'k')
        hold on
        plot(hs2,m(3,:),alpha,'ob')
        plot(hs2,[0 261770],[mean(alpha)-messages.setup.wavedet.thfraction mean(alpha)-messages.setup.wavedet.thfraction],':g')
        plot(hs2,[0 261770],[th th],':r'),
        plot(hs2,m(3,ind),alpha(ind),'xr')
        set(hs2,'xlim',aa)
        ylabel(hs2,'\alpha')
        
        %legend('alpha','alpha','mean alpha-thfraction','th','Location','NorthEastOutside')
        
    end
    if ~isempty(ind)			     % and high frequency noise peaks
        m(:,ind) = [];
        signo(ind)=[];
    end
    thresinterval = ceil(messages.setup.wavedet.intvlthr2* messages.setup.wavedet.freq);
    % threshold interval to consider isolated maximum lines
    ind = find( ((m(1,2:end-1)-m(1,1:end-2))>thresinterval) ...
        &      (m(1,3:end)- m(1,2:end-1))>thresinterval)+1;
    
    if size(m,2)<2 && ~isempty(ind), ind=1; end  % when only one maximum %16DEZ08
    if ff==1
%         axes(fig1);
        plot(figax1,m(1,ind),w(m(1,ind)),'+r')
        if ~isempty(ind)
            ll=[ll {'isolated mm lines'}];
        end
%         axes(fig1);
        legend(figax1,['WT_1',ll],'Location','NorthEastOutside')
%         axes(fig3);
        plot(figax3,m(3,ind),w(m(3,ind),3),'+r')
        legend(['WT_3',ll],'Location','NorthEastOutside')
        
        try
%         axes(figECGWT1);
        plot(figWT1ECG,m(1,ind),w(m(1,ind),1),'+r')
        plot(figWT1ECG,[m(1,ind(1))-thresinterval m(1,ind(1))-thresinterval],[-5000 5000],':r')
        plot(figWT1ECG,[m(1,ind(1))+thresinterval m(1,ind(1))+thresinterval],[-5000 5000],':r')
%         axes(figECGWT2);
        plot(figWT2ECG,m(2,ind),w(m(2,ind),2),'+r')
%         axes(figECGWT3);
        plot(figWT3ECG,m(3,ind),w(m(3,ind),3),'+r')
        catch me
            me.message = 'Error plotting WT signal';
        end
    end
    m(:,ind) = [];                     %Discard isolated maximum lines
    signo(ind)=[];
    % Threshold interval for considering redundancy
    redundant= [];
    thresinterval = ceil(messages.setup.wavedet.intvlthr1*messages.setup.wavedet.freq);           % 120 ms. (li)
    for l = find(signo>0),            % For each positive maximum line
        if ~any(redundant ==l),         % If it has not been declared redundant yet
            ind=find((m(3,:)>m(3,l)-messages.setup.wavedet.intvlthr1_2*thresinterval)&(m(3,:)<m(3,l)+messages.setup.wavedet.intvlthr1_2*thresinterval)&signo>0);
            % index of positive lines near the present one (including it)
            if length(ind)>1,              % If more than one --> redundancy
                [mx,ind2]= max(abs(w(m(3,ind),3))); %#ok<ASGLU>
                ind(ind2)=[];
                redundant = [redundant ind];   %#ok<AGROW> % All but the greatest are redundant
            end
        end
    end
    if ff==1
%         axes(fig1);
        plot(figax1,m(1,redundant),w(m(1,redundant)),'sm')
%         axes(fig3);
        
        plot(figax3,m(3,redundant),w(m(3,redundant),3),'sm')
        if ~isempty(redundant)
            ll=[ll {'redundant mm lines'}];
        end
%         axes(figECGWT1);
        plot(figWT1ECG,m(1,redundant),w(m(1,redundant)),'sm')
%         axes(figECGWT2);
        plot(figWT2ECG,m(2,redundant),w(m(2,redundant),2),'sm')
%         axes(figECGWT3);
        plot(figWT3ECG,m(3,redundant),w(m(3,redundant),3),'sm')
        plot(figWT3ECG,[m(3,206)-messages.setup.wavedet.intvlthr1_2*thresinterval m(3,206)-messages.setup.wavedet.intvlthr1_2*thresinterval],[-5000 5000],':r')
        plot(figWT3ECG,[m(3,206)+messages.setup.wavedet.intvlthr1_2*thresinterval m(3,206)+messages.setup.wavedet.intvlthr1_2*thresinterval],[-5000 5000],':r')

    end
    m(:,redundant)=[];		  % Discard redundant lines
    signo(redundant)=[];
    redundant =[];
    for l = find(signo>0),		  % For each remaining positive maximum line
        ind=find((m(3,:)>m(3,l)-thresinterval)&(m(3,:)<m(3,l)+thresinterval)&signo<0);
        % Search for negative minima near it
        if length(ind)>1,             % If more than one ---> redundancy
            aux = abs(w(m(3,ind),3)./(m(3,ind)'-m(3,l)));
            % auxiliary variable: height over distance
            [mx,indmx]=max(aux);
            ind2=ind;
            aux(indmx)=[]; ind2(indmx)=[];
            if all((mx./aux)>peakcriteria),      % RULE 2 (see PFC or Li's paper)
                redundant = [redundant ind2]; %#ok<AGROW>
            else
                [aux,aux2]=min(abs(m(3,l)-m(3,ind))); %#ok<ASGLU>
                ind(aux2)=[];		   % RULE 1 (see PFC or Li's paper)
                redundant = [redundant ind]; %#ok<AGROW>
            end
        end
    end
    if ff==1
%         axes(fig1);
        plot(figax1,m(1,redundant),w(m(1,redundant)),'*r')
%         axes(fig3);
        plot(figax3,m(3,redundant),w(m(3,redundant),3),'*r')
        if ~isempty(redundant)
            ll=[ll {'redundant mm lines'}];
        end
        try
%             axes(figECGWT1);
            plot(figWT1ECG,m(1,redundant(end)),w(m(1,redundant(end))),'sm')
%             axes(figECGWT2);
            plot(figWT2ECG,m(2,redundant(end)),w(m(2,redundant(end)),2),'sm')
%             axes(figECGWT3);
            plot(figWT3ECG,m(3,redundant(end)),w(m(3,redundant),3),'sm')
            plot(figWT3ECG,[m(3,182)-messages.setup.wavedet.intvlthr1_2*thresinterval m(3,182)-messages.setup.wavedet.intvlthr1_2*thresinterval],[-5000 5000],':r')
            plot(figWT3ECG,[m(3,182)+messages.setup.wavedet.intvlthr1_2*thresinterval m(3,182)+messages.setup.wavedet.intvlthr1_2*thresinterval],[-5000 5000],':r')
        catch me
            me.message = 'Error plotting WT signal';
        end
    end
    m(:,redundant)=[];		  % Discard redundant lines
    signo(redundant)=[];
    redundant =[];
    for l = find(signo<0),            % For each remaining negative minimum line
        ind = find((m(3,:)>m(3,l)-thresinterval)&(m(3,:)<m(3,l)+thresinterval)&signo>0);
        % Search for positive maxima near it
        if length(ind)>1,		   % If more than one ----> redundancy
            aux = abs(w(m(3,ind),3)./(m(3,ind)'-m(3,l)));
            % auxiliary variable: height over distance
            [mx,indmx]=max(aux);
            ind2=ind;
            aux(indmx)=[]; ind2(indmx)=[];
            if all((mx./aux)>peakcriteria),
                redundant = [redundant ind2]; %#ok<AGROW> % RULE 2
            else
                [aux,aux2]=min(abs(m(3,l)-m(3,ind))); %#ok<ASGLU>
                ind(aux2)=[];		        % RULE 1
                redundant = [redundant ind]; %#ok<AGROW>
            end
        end
    end
    if ff==1
%        axes(fig1);
        plot(figax1,m(1,redundant),w(m(1,redundant)),'*r')
%        axes(fig3);
        plot(figax3,m(3,redundant),w(m(3,redundant),3),'*r')
        if ~isempty(redundant)
            ll=[ll {'redundant mm lines'}];
        end
    end
    m(:,redundant)=[];	            % Discard redundant lines
    signo(redundant)=[];
    %%%% isolated maximum lines resulting from Discarded redundant
    ind = find( ((m(1,2:end-1)-m(1,1:end-2))>thresinterval) ...
        &      (m(1,3:end)- m(1,2:end-1))>thresinterval)+1;
    if size(m,2)<2 && ~isempty(ind), ind=1; end  % when only one maximum %16DEZ08
    if ff==1
%        axes(fig1);
        legend(figax1,['WT_1',ll],'Location','NorthEastOutside')
%        axes(fig3);
        legend(figax3,['WT_3',ll],'Location','NorthEastOutside')
%        axes(fig1);
        plot(figax1,m(1,ind),w(m(1,ind)),'+m')
%        axes(fig3);
        plot(figax3,m(3,ind),w(m(3,ind),3),'+m')
        if ~isempty(ind)
            ll=[ll {'isolated mm lines'}];
        end
%        axes(fig1);
        legend(figax1,['WT_1',ll],'Location','NorthEastOutside')
%        axes(fig3);
        legend(figax3,['WT_3',ll],'Location','NorthEastOutside')
    end
    m(:,ind) = [];                     %Discard isolated maximum lines
    signo(ind)=[];
    if size(m,2)<2
        m=[];
        signo=[];
    end
    eliminar=[];
    %%%%extra protection% 23MAR09
    for ii=1:size(m,2)
        pa = picant(w(max(1,m(2,ii)-round(messages.setup.wavedet.pictime*messages.setup.wavedet.freq)):m(2,ii),2),m(2,ii));
        % first peak before detected qrs position at scale 2
        pp = picpost(w(m(2,ii):min(size(w,1),m(2,ii)+round(messages.setup.wavedet.pictime*messages.setup.wavedet.freq)),2),m(2,ii));
        %if isempty(pa) || isempty(pp)
        if isempty(pa) && isempty(pp) %NOV2011
            eliminar=[eliminar ii]; %#ok<AGROW>
        end
    end
    if ff==1
%        axes(fig1);
        plot(figax1,m(1,eliminar),w(m(1,eliminar)),'.r')
%        axes(fig3);
        plot(figax3,m(3,eliminar),w(m(3,eliminar),3),'.r')
        if ~isempty(eliminar)
            ll=[ll {'extra protection'}];
        end
%        axes(fig1);
        legend(figax1,['WT_1',ll],'Location','NorthEastOutside')
%        axes(fig3);
        legend(figax3,['WT_3',ll],'Location','NorthEastOutside')
    end
    m(:,eliminar)=[];
    signo(eliminar)=[];
     % QRS peak detection / wavelet Zero cross detection
    timelap = ceil(messages.setup.wavedet.timelapthr*messages.setup.wavedet.freq);   % the two maximum lines defining a QRS
    % should not be separated more than timelapthr
    time = [];
    aux=[];
    auxm=[];
    if size(m,2)>1 % 18 OUT2011
        for l = find(signo>0), % For each positive maximum line
            if (l==1)                         % Special case: first line
                if (signo(2)<0)&&(m(1,2)-m(1,1)<thresinterval),
                    ind = zerocros(w(m(1,l):min(m(1,l)+timelap,size(w,1)),1));
                    % Zero crossing at scale 1
                    time = [time ind+m(1,l)-1]; %#ok<AGROW>
                    aux=[aux abs(w(m(1,l))-w(m(1,l+1)))]; %#ok<AGROW>
                    auxm=[auxm; [m(1,l) m(1,l+1)]]; %#ok<AGROW>
                end
            elseif (l==size(m,2))             % Special case: last line
                if (signo(l-1)<0)&&(m(1,end)-m(1,end-1)<thresinterval),
                    ind = zerocros(w(m(1,l-1):min(m(1,l-1)+timelap,size(w,1)),1));
                    time = [time ind+m(1,l-1)-1]; %#ok<AGROW>
                    aux=[aux abs(w(m(1,l-1))-w(m(1,l)))]; %#ok<AGROW>
                    auxm=[auxm ; [m(1,l-1) m(1,l)]]; %#ok<AGROW>
                end
            elseif signo(l+1)<0 && ((signo(l-1)>0) ...
                    || ((m(1,l+1)-m(1,l))<(m(1,l)-m(1,l-1))))
                ind = zerocros(w(m(1,l):min(m(1,l)+timelap,size(w,1)),1));
                time = [time ind+m(1,l)-1]; %#ok<AGROW>
                aux=[aux abs(w(m(1,l))-w(m(1,l+1)))]; %#ok<AGROW>
                auxm=[auxm;[ m(1,l) m(1,l+1)]]; %#ok<AGROW>
            elseif signo(l-1)<0, 
                ind = zerocros(w(m(1,l-1):min(m(1,l-1)+timelap,size(w,1)),1));
                time = [time ind+m(1,l-1)-1]; %#ok<AGROW>
                aux=[aux abs(w(m(1,l-1))-w(m(1,l)))]; %#ok<AGROW>
                auxm=[auxm;[ m(1,l-1) m(1,l)]]; %#ok<AGROW>
            end
        end
    end
    if ff==1
       figax2 = axes(fig2);
        plot(figax2,time,messages.sig(time),'ob')
        ll_ECG={'ECG','QRS candidates'};
        legend(figax2,ll_ECG,'Location','NorthEastOutside')
    end
    % Refractary period after a QRS detection
    if ~isempty(time),
        rr = (time(2:end)-time(1:end-1));
        ind = find(rr<ceil(refrper*messages.setup.wavedet.freq));             
         for auxi=1:length(ind) %RUTE 27Jun11
             [M,ii]=min([aux(ind(auxi)) aux(ind(auxi)+1)]); %#ok<ASGLU> % 18JUL2011
             if ii==2, ind(auxi)=ind(auxi)+1; end
         end
%       ind=ind+1; % olde version % always the second one is eliminated!!!
        if ff==1
%            axes(fig2);
            plot(figax2,time(ind),messages.sig(time(ind)),'xr')
%            axes(fig1);
            plot(figax1,m(1,ind),w(m(1,ind)),'xr')
%            axes(fig3);
            plot(figax3,m(3,ind),w(m(3,ind),3),'xr')
            if ~isempty(ind)
                ll=[ll {'rejected by the refractary period'}];
                ll_ECG=[ll_ECG {'rejected by the refractary period'}];
            end
%            axes(fig2);
            legend(figax2,ll_ECG,'Location','NorthEastOutside')
%            axes(fig1);
            legend(figax1,['WT_1',ll],'Location','NorthEastOutside')
%            axes(fig3);
            legend(figax3,['WT_3',ll],'Location','NorthEastOutside')
        end
        time(ind)=[]; %RUTE 27Jun11
        aux(ind)=[];
        if abs(samp(time(1)) - lastqrs)< refrper*messages.setup.wavedet.freq,   % Refractory period for
            time(1)=[];
            aux(1)=[];%#ok<NASGU> % last beat of previous
        end        % fragment       
    end
    % search back if no detection in messages.setup.wavedet.rrmaxlim RR median
    if ~isempty(time),
        %%%%%%%%%%%%%%%% caso em que rr tem dimensao inferior ao necessario
        if length(time)>1  %%%%%%%%%%%%%%%%introduzido a 17/05/02 Rute
            rr = ([max(samp(time(1))-lastqrs,time(1)) time(2:end)-time(1:end-1)]); %rr = ([samp(time(1))-lastqrs, time(2:end)-time(1:end-1)]);
            rrmed = [rr(1) 0.5*(rr(1)+rr(2)) median([rr(3:end);rr(2:end-1);rr(1:end-2)])];
            ind = find((rr(2:end)>messages.setup.wavedet.rrmaxlim*rrmed(1:end-1)) & (rr(1:end-1)>messages.setup.wavedet.rrmaxmaxlim*messages.setup.wavedet.freq) |  rr(2:end)>messages.setup.wavedet.rrmaxlimsec*messages.setup.wavedet.freq);
            %%%% 22NOV2011
            if messages.setup.wavedet.refrper_reduc>0 && refrper>0.250 % 22NOV2011
                while sum(rrmed<refrper*messages.setup.wavedet.freq*4/3)/length(rrmed)>0.6
                    refrper=refrper*messages.setup.wavedet.refrper_reduc;
                end % 22NOV2011
                if messages.setup.wavedet.refrper(end)~=refrper
                    messages.setup.wavedet.refrper(end+1)=refrper;
                    messages.setup.wavedet.refrperdef(end+1)=samp(time(1));
                end
            end
            if ff==1
                if ~isempty(ind)
                    ll=[ll {'searchback'}];
                    ll_ECG=[ll_ECG {'searchback'}];
%                    axes(fig2);
                    plot(figax2,time(ind),messages.sig(time(ind))+50,'vc')
                    legend(ll_ECG,'Location','NorthEastOutside')
%                    axes(fig1);
                    plot(figax1,time(ind),w(time(ind),1)+200,'vc')
%                    axes(fig3);
                    plot(figax3,time(ind),w(time(ind),3)+200,'vc')
%                    axes(fig1);
                    legend(figax1,'WT_1',ll,'Location','NorthEastOutside')
%                    axes(fig3);
                    legend(figax3,'WT_3',ll,'Location','NorthEastOutside')
                end
            end
            for l = ind,               % For all rr intervals greater than 1.5RR
                interv = (time(l)+ceil(refrper*messages.setup.wavedet.freq)) : (time(l+1)-ceil(refrper*messages.setup.wavedet.freq));
                if length(interv)>messages.setup.wavedet.rrmax_rbt*messages.setup.wavedet.freq,  % for robustness
                   % [nuevo, messages] = searchbk(w(interv,:),interv(1),eps(:,nlead)',messages.setup.wavedet.freq,messages.setup.wavedet.th_schb,messages);
                   [nuevo, messages] = searchbk(w(interv,:),interv(1),eps(:,nlead)',messages.setup.wavedet.freq,th,messages); %FEB16_2012
                    % searchbk returns the new time of occurrence of a QRS
                    time = [time nuevo];     %#ok<AGROW> % New times are added
                    %  time = sort(time);    	     % and sorted
                    if ff==1
%                        axes(fig2);
                        plot(figax2,nuevo,messages.sig(nuevo),'oc')
                        if ~isempty(nuevo)
                            ll_ECG=[ll_ECG {'new candidates by searchback'}]; %#ok<AGROW>
                        end
                        legend(figax2,ll_ECG,'Location','NorthEastOutside')
                    end
                end
            end %%%%%%%%%%%%%%%%introduzido a 17/05/02 Rute
        else % only one detection DEZ2011
            if time > max(messages.setup.wavedet.gapthrs,ceil(refrper*messages.setup.wavedet.freq))%initial gap % 24ABRIL2010
%             if time >  messages.setup.wavedet.gapthrs %initial gap
                if (~exist('th','var') || isnan(th)), th = messages.setup.wavedet.thmax ; end
                [nuevo, messages]= searchbk(w(1:(time-ceil(refrper*messages.setup.wavedet.freq))),1,eps(:,nlead)',messages.setup.wavedet.freq,th,messages);
            time = [time nuevo]; 
            end
        end
    else  % in case there is no detection, searchback anyway
        if (~exist('th','var') || isnan(th)), th = messages.setup.wavedet.thmax ; end
        %nuevo = searchbk(w,1,eps,messages.setup.wavedet.freq, messages.setup.wavedet.th_schb); %Rute 02.Dec.04
        [nuevo, messages]= searchbk(w,1,eps(:,nlead)',messages.setup.wavedet.freq, th,messages);
        time = nuevo;
        if ff==1
%            axes(fig2);
            plot(figax2,nuevo,messages.sig(nuevo),'oc')
            if ~isempty(nuevo)
                ll_ECG=[ll_ECG {'new candidates by searchback'}];
            end
            legend(figax2,ll_ECG,'Location','NorthEastOutside')
        end
    end
    if ~isempty(time),
        time = sort(time);  % Rute 14DEz
        if length(time)>1
            rr = ([samp(time(1))-lastqrs, time(2:end)-time(1:end-1)]);
            rrmed = [rr(1) 0.5*(rr(1)+rr(2)) median([rr(3:end);rr(2:end-1);rr(1:end-2)])]; %RUTE OUT2011
        else
            rrmed =1;
        end
        if length(w)-time(end) >  messages.setup.wavedet.gapthrs*rrmed,  % If there is a final "gap"
            interv= (time(end)+ceil(refrper*messages.setup.wavedet.freq) : length(w)); %RUTE OUT2011
            if length(interv)>messages.setup.wavedet.rrmax_rbt*messages.setup.wavedet.freq,
                [nuevo, messages] = searchbk(w(interv,:),interv(1),eps(:,nlead)',messages.setup.wavedet.freq, th,messages);
                time = [time nuevo];
                if ff==1
%                    axes(fig2);
                    plot(figax2,nuevo,messages.sig(nuevo),'oc')
                    if ~isempty(nuevo)
                        ll_ECG=[ll_ECG {'new candidates by searchback'}];
                    end
                    legend(figax2,ll_ECG,'Location','NorthEastOutside')
                end
            end
        end
        %%% The first and last beats must be complete, to detect P-wave and T-wave %%%
        %%% We use an overlap that assures that the beat will be complete always in
        %%% one of the segments.
        if (time(1) <= messages.setup.wavedet.pwaveout*messages.setup.wavedet.freq),
            time(1)=[];
        end
        % If the P-wave is not in the present segment, discard the beat, because the
        % the beat should have been detected in the last segment.
        if ~isempty(lastqrs)
            time(samp(time)<= lastqrs + refrper*messages.setup.wavedet.freq)=[];
        end
        if ~isempty(time) %%%%%%% Rute 06/05/02 caso time(indrep)=time
            % if beat before ultimo_anot in the previous segment %6AGO2007
            if samp(time(1))<=ultimo_anot
                time(1)=[];
            end           
            if ~isempty(time) && (length(w)-time(end) < messages.setup.wavedet.twaveout*messages.setup.wavedet.freq),  % si no ye a onda T drento
                time(end)=[];                               % de sig
            end
            % If last beat's T-wave is not in 'sig', the beat is detected in next one
            if ~isempty(time) %Rute 10.04.05
                timeqrs = [timeqrs samp(time)];  % QRS times
                lastqrs = samp(time(end));
            end
        end
        intervalo = [numlatdet+1 numlatdet+length(time)];
        % numeration of the beats processed in this segment
        numlatdet = intervalo(2);
        position.qrs(intervalo(1):intervalo(2)) = samp(time);  % Fill QRS position
    end
    if ff==1
%        axes(fig2);
        plot(figax2,time,messages.sig(time),'.g')
        ll_ECG=[ll_ECG {'final SL marks'}];
        legend(figax2,ll_ECG,'Location','NorthEastOutside')
%        axes(fig1);
        legend(figax1,['WT_1',ll],'Location','NorthEastOutside')
%        axes(fig3);
        legend(figax3,['WT_3',ll],'Location','NorthEastOutside')
    end    
else
    intervalo=[];time=[]; %Rute 13Set06
end

