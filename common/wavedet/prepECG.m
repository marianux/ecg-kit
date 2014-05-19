function [error_code error_mess sig_ECGout signal_levels sig_level baseline_level ac_level baseline_correlation]=prepECG(sig_ECGin);


 %-------------------------------------------
% preFILT1 - BASELINE AND AC FILTER DATA (ENTIRE SIGNAL)
%-------------------------------------------

sig_ECGout=zeros(size(sig_ECGin));
bfiltfrek=0.9;  % maybe higher for short signals
[B1,A1]=butter(3,bfiltfrek/50);
[B2,A2]=butter(4,[48 52]/500,'stop');% maybe add harmonics
[B3,A3]=butter(4,[58 62]/500,'stop');% maybe add harmonics

block_AMP=2000;
block_BASE=100;

dylenAMP = floor(size(sig_ECGout,2)/block_AMP);
dylenBASE = floor(size(sig_ECGout,2)/block_BASE);

Ecg100 = zeros(size(sig_ECGin,1), ceil(length(sig_ECGin(1,1:2:end))/5));
baseline_level = zeros(1,size(sig_ECGin,1));
ac_level = zeros(size(sig_ECGin,1),1);
sig_level = zeros(size(sig_ECGin,1),4);
leads=size(sig_ECGin,1);

for lead=1:size(sig_ECGin,1)
    Ecg100(lead,:)=decimate(sig_ECGin(lead,1:2:end),5);
    Baseline=filtfilt(B1,A1,Ecg100(lead,:));
    Baseline=interp(Baseline,10);
    baseline_level(lead)=std(Baseline);
    sig_ECGout(lead,:)=sig_ECGin(lead,:)-Baseline(1:size(sig_ECGin,2));
    tmp=sig_ECGout(lead,:);
    sig_ECGout(lead,:)=filtfilt(B2,A2,sig_ECGout(lead,:));       %
    ac_level(lead,1)=std(sig_ECGout(lead,:)-tmp);
    tmp=sig_ECGout(lead,:);
    sig_ECGout(lead,:)=filtfilt(B3,A3,sig_ECGout(lead,:));       %
    ac_level(lead,2)=std(sig_ECGout(lead,:)-tmp);

    DynamicsP=zeros(1,dylenAMP);
    DynamicsN=zeros(1,dylenAMP);
    for i=1:floor(size(sig_ECGout,2)/block_AMP)
        idx = (i-1)*block_AMP+1:i*block_AMP;
        DynamicsP(i)=max(sig_ECGout(lead,idx));
        DynamicsN(i)=min(sig_ECGout(lead,idx));         
    end;
    
    DynbaseP=zeros(1,dylenBASE);
    DynbaseN=zeros(1,dylenBASE);
    for i=1:floor(size(sig_ECGout,2)/block_BASE)
        idx = (i-1)*block_BASE+1:i*block_BASE;
        DynbaseP(i)=max(sig_ECGout(lead,idx));
        DynbaseN(i)=min(sig_ECGout(lead,idx));        
    end;
    Dynbase=DynbaseP-DynbaseN;
    
    tmpz=zeros(size(sig_ECGout(lead,:)));
           tmpz(sig_ECGout(lead,:)>3*median(DynamicsP(find(DynamicsP>1))))=1;
           tmpz(sig_ECGout(lead,:)<3*median(DynamicsN(find(abs(DynamicsN)>1))))=1;
           tmpz=filtfilt(1/175.* ones(1,175),1,tmpz);
          sig_ECGout(lead,tmpz>0)=0;
    
    
    
    sig_level(lead,1)=median(DynamicsP(find(DynamicsP>1)));
    sig_level(lead,2)=median(DynamicsN(find(abs(DynamicsN)>1)));

    
    if sum(DynbaseP>1)>0
        sig_level(lead,3)=prctile(DynbaseP(find(DynbaseP)),20);
    else
       sig_level(lead,3)=NaN;
    end;   
   if sum(Dynbase>1)>0    
    sig_level(lead,4)=prctile(Dynbase(find(Dynbase>1)),20);
    else
       sig_level(lead,4)=NaN;
    end;      
    
    sig_ECGout(lead,:)=filtfilt(0.05*ones(20,1),1,sig_ECGout(lead,:));
    
    if lead==1
        sig_ECGcut=sig_ECGout(1,:);
        sig_ECGcut(sig_ECGcut>sig_level(1,3))=sig_level(1,3);
        sig_ECGcut(sig_ECGcut<(sig_level(1,3)-sig_level(1,4)))=sig_level(1,3)-sig_level(1,4);

        Y1=xcorr(sig_ECGcut,'unbiased');
        [MaxPos1 MaxVal1]=locmax(Y1);
        MaxPos1=MaxPos1(MaxVal1>0);
        MaxVal1=MaxVal1(MaxVal1>0);

         

        if length(MaxPos1)>1
            baseline_correlation.maincycle=MaxPos1(ceil(length(MaxPos1)/2)+1)-MaxPos1(ceil(length(MaxPos1)/2));%Avstånd mellan main peak och second peak i corrsig.
            baseline_correlation.maincyclecorr=MaxVal1(ceil(length(MaxPos1)/2)+1)/MaxVal1(ceil(length(MaxPos1)/2)); %Amplitud ratio
            baseline_correlation.partclosetomc=length(find(abs(diff(MaxPos1)-baseline_correlation.maincycle)<20))/length(diff(MaxPos1)); % andel av korrpeakar som passar med flength.
        else
            baseline_correlation.maincycle=NaN;
            baseline_correlation.maincyclecorr=NaN;
            baseline_correlation.partclosetomc=NaN;
        end;

        
        
    end; %lead==1
end;%lead

signal_levels=[mean(baseline_level) mean(ac_level,1) mean(sig_level,1) sig_level(1,3) sig_level(1,4)];

leads_normalbeatrange=find(sig_level(1:min(leads,3),1)<10000 & sig_level(1:min(leads,3),2)>-10000 & (sig_level(1:min(leads,3),1)-sig_level(1:min(leads,3),2))<12000 & (sig_level(1:min(leads,3),1)-sig_level(1:min(leads,3),2))>200);
if length(leads_normalbeatrange)<min(leads,3);
    error_code=103;
    error_mess='103 (A.2) - Unexpected beat range (some of leads 1-3)';
    return
        end;
if (sig_level(1:min(leads,3),3)-sig_level(1:min(leads,3),4))>0.5*(sig_level(1:min(leads,3),1)-sig_level(1:min(leads,3),2))
    error_code=104;
    error_mess='104 (A.2) - Baseline larger than 50% of beat range (some of leads 1-3)';
    return
end;
error_code=0;
error_mess='';   
  

