function [error_code error_mess sig_ECG,lab_ECG,Leads,samples,e_ECG]=readECG_v1(sampling_rate,data,timelim)
e_ECG = [];
sig_ECG = double(data);
%sig_ECG = sig_ECG(2:end);
lab_ECG = '';

% Error 101
if isempty(sig_ECG)
    Leads=0;
    TotLength=0;  
    error_code=101; 
    error_mess='101 (A.1) - Empty signal';    
    return
end;
 
 
 % Sampla om om det behövs
if sampling_rate~=1000;
    if sampling_rate==125;
        up=8;
        down=1;
    elseif sampling_rate==128;
        up=125;
        down=16;
    elseif sampling_rate==256;    
        up=125;
        down=32;
    elseif sampling_rate==250;
        up=4;
        down=1;
    elseif sampling_rate==512    
        up=125;
        down=64;
    elseif sampling_rate==500    
        up=2;
        down=1;
    elseif sampling_rate==1024    
        up=125;
        down=128;  
    elseif sampling_rate==2000    
        up=1;
        down=2;  
    else
        up=1000;
        down=sampling_rate;
    end;
else
    up=1;
    down=1;
end;    


if up~=1 || down~=1
    
% changed by MEC
%     sig_ECGtmp=zeros(size(sig_ECG,1),ceil(size(sig_ECG,2)*up/down));
%     for j=1:size(sig_ECG,1)
%         sig_ECGtmp(j,:)=decimate(interp(sig_ECG(j,:),up),down); 
%     end; 
%     sig_ECG=sig_ECGtmp;
%     clear sig_ECGtmp;
    [P Q] = rat(up/ down);
    sig_ECG = resample(sig_ECG', P, Q)';
end;

% Error 102
if size(sig_ECG,2)<10000
    Leads=0;
    TotLength=0; 
    error_code=102; 
    error_mess='102 (A.1) - Signal<10sec';    
    return
end;
 
% Info about read signal
[Leads, samples]=size(sig_ECG);

% Temporary length limitation
sig_ECG=sig_ECG(:,1:min(timelim,samples));
samples=size(sig_ECG,2); %if it was cut


 % Spara invasiva mätningar separat
 if Leads>12
    e_ECG=sig_ECG(13:end,:);
 end;    
 
 % Spara bara EKG is sig_ECG
 Leads=min(Leads,12);
 sig_ECG=sig_ECG(1:Leads,:);



 
% Gallring klar   
switch Leads
    case 12
        lab_ECG=['V1 ';'V2 ';'V3 ';'V4 ';'V5 ';'V6 ';'I  ';'II ';'III';'aVR';'aVL';'aVF'];
    case 9
        lab_ECG=['V1 ';'V2 ';'V3 ';'V4 ';'V5 ';'V6 ';'I  ';'II ';'III'];
    case 8
        lab_ECG=['V1 ';'V2 ';'V3 ';'V4 ';'V5 ';'V6 ';'I  ';'II '];
    case 6
        lab_ECG=['V1 ';'V2 ';'V3 ';'V4 ';'V5 ';'V6 '];
    case 3
        lab_ECG=['1/X';'2/Y';'3/Z'];
    case 2
        lab_ECG=['1  ';'2  '];        
    case 1
        lab_ECG='1  ';
end

lab_ECG(lab_ECG==' ')='_'; % space does not work in variables



error_code=0;
error_mess='';    

