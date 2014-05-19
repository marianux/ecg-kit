function [timeqrs,indexes,tol,messages]=qrscandidatesnew(timeqrs1,timeqrs2,timeqrs3,heasig,tol,messages)
% [timeqrs,indexes]=qrscandidates(timeqrs1,timeqrs2,timeqrs3,heasig,tol)
%
% Construction QRS positions candidates vector from 3 single detections vectors
% a QRS is admitted as candidate if:
%   - it was found in 3 leads and the middle mark differ less than tol from the others
%   - it was found in 2 leads and the marks differ less than tol 
%
% marks that not fullfill the above are replaced by NaN and a QRS positions
% candidates vector is constructed with the detections in each lead in one
% line, one QRS in each column
%
% INPUT:
% timeqrs1,timeqrs2,timeqrs3  -  line vector swith single lead detections
%                                (in samples)
% heasig - header information
% tol - maximum distance to be admited as same complex (in sec)
%       by default half of the default refractary period tol = 0.275/2 sec
% 
% OUTPUT:
% timeqrs - QRS positions candidates vector
% indexes - beats corresponding to the QRS positions candidates vector
%
% Rute Almeida  14.APR.2005
% Last update: 26JUL2011
%
% MATLAB Version R13
if nargin<6
    messages.status=1;
end
if ~isfield(messages,'warnings'), messages.warnings=[]; end
spf=heasig.spf_ecg; %24.MAR.09
if nargin==4
    tol=ceil(0.275/2*heasig.freq*spf); %24.MAR.09
elseif nargin<4
    messages.errors=[messages.errors {'Fatal error in qrscandidatesnew: not enough inputs.'}];
    warning(char(messages.errors(end)))
    messages.errors_desc=[messages.errors_desc 'Mandatory inputs not defined.'];
    messages.status=0;
    return
end
indexes=[];
indexes1=1:length(timeqrs1);
indexes2=1:length(timeqrs2);
indexes3=1:length(timeqrs3);
timeqrs=[];
lengthaux=([length(timeqrs1) length(timeqrs2) length(timeqrs3)]);

for g=1:min(lengthaux)
    [aux auxi]=sort([timeqrs1(g) timeqrs2(g) timeqrs3(g)]);
    if aux(1)>(aux(2)-tol-1)
        if (aux(2)+tol+1)>aux(3)
            timeqrs=[timeqrs [timeqrs1(g);timeqrs2(g);timeqrs3(g)]]; %#ok<AGROW>
            indexes=[indexes [indexes1(g);indexes2(g);indexes3(g)]]; %#ok<AGROW>
        else
            if auxi(3)==1
                timeqrs1(g+1:end+1)=timeqrs1(g:end);
                indexes1(g+1:end+1)= indexes1(g:end);
                indexes1(g)=NaN; 
                timeqrs1(g)=NaN;
            else if auxi(3)==2
                    timeqrs2(g+1:end+1)=timeqrs2(g:end);
                    indexes2(g+1:end+1)= indexes2(g:end);
                    timeqrs2(g)=NaN;
                    indexes2(g)=NaN;
                else
                    timeqrs3(g+1:end+1)=timeqrs3(g:end);
                    indexes3(g+1:end+1)= indexes3(g:end);
                    timeqrs3(g)=NaN;
                    indexes3(g)=NaN;
                end
            end
            timeqrs=[timeqrs [timeqrs1(g);timeqrs2(g);timeqrs3(g)]]; %#ok<AGROW>
            indexes=[indexes [indexes1(g);indexes2(g);indexes3(g)]]; %#ok<AGROW>
        end
    else
        
        if auxi(1)==1
            g1=g+1;
            flag=0;
            while  g1<=length(timeqrs1) && flag==0 && timeqrs1(g1)<(aux(2)+tol+1)
                if (timeqrs1(g1)>(aux(2)-tol-1))
                    flag=1;
                    timeqrs1(g:end)=[timeqrs1(g1:end) NaN*ones(1,(g1-g))];
                    indexes1(g:end)=[indexes1(g1:end) NaN*ones(1,(g1-g))];
                end
                g1=g1+1;
            end
            if flag==0       
                if g1>length(timeqrs1)
                    timeqrs1(g:end)= NaN*ones(1,length(timeqrs1(g:end)));
                    indexes1(g:end)= NaN*ones(1,length(timeqrs1(g:end)));
                else
                    timeqrs1(g:end)=[timeqrs1((g1-1):end) NaN*ones(1,(g1-1-g))];
                    indexes1(g:end)=[indexes1((g1-1):end) NaN*ones(1,(g1-1-g))];
                    timeqrs1(g)=NaN;
                    indexes1(g)=NaN;
                end
            end
        elseif auxi(1)==2
            g1=g+1;
            flag=0;
            while  g1<=length(timeqrs2) && flag==0 && timeqrs2(g1)<(aux(2)+tol+1)
                if (timeqrs2(g1)>(aux(2)-tol-1))
                    flag=1;
                    timeqrs2(g:end)=[timeqrs2(g1:end) NaN*ones(1,(g1-g))];
                    indexes2(g:end)=[indexes2(g1:end) NaN*ones(1,(g1-g))];
                end
                g1=g1+1;
            end
            if flag==0       
                if g1>length(timeqrs1)
                    timeqrs2(g:end)= NaN*ones(1,length(timeqrs2(g:end)));
                    indexes2(g:end)= NaN*ones(1,length(timeqrs2(g:end)));
                else
                    timeqrs2(g:end)=[timeqrs2((g1-1):end) NaN*ones(1,(g1-1-g))];
                    indexes2(g:end)=[indexes2((g1-1):end) NaN*ones(1,(g1-1-g))];
                    timeqrs2(g)=NaN;
                    indexes2(g)=NaN;
                end
            end          
        else
            g1=g+1;
            flag=0;
            while  g1<=length(timeqrs3) && flag==0 && timeqrs3(g1)<(aux(2)+tol+1)
                if (timeqrs3(g1)>(aux(2)-tol-1))
                    flag=1;
                    timeqrs3(g:end)=[timeqrs3(g1:end) NaN*ones(1,(g1-g))];
                    indexes3(g:end)=[indexes3(g1:end) NaN*ones(1,(g1-g))];
                end
                g1=g1+1;
            end
            if flag==0       
                if g1>length(timeqrs1)
                    timeqrs3(g:end)= NaN*ones(1,length(timeqrs3(g:end)));
                    indexes3(g:end)= NaN*ones(1,length(timeqrs3(g:end)));
                else
                    timeqrs3(g:end)=[timeqrs3((g1-1):end) NaN*ones(1,(g1-1-g))];
                    indexes3(g:end)=[indexes3((g1-1):end) NaN*ones(1,(g1-1-g))];
                    timeqrs3(g)=NaN;
                    indexes3(g)=NaN;
                end
            end
        end
        
        if aux(3)-aux(2)>(tol+1)
            if auxi(3)==1
                timeqrs1(g+1:end+1)=timeqrs1(g:end);
                indexes1(g+1:end+1)= indexes1(g:end);
                indexes1(g)=NaN; 
                timeqrs1(g)=NaN;
            else if auxi(3)==2
                    timeqrs2(g+1:end+1)=timeqrs2(g:end);
                    indexes2(g+1:end+1)= indexes2(g:end);
                    timeqrs2(g)=NaN;
                    indexes2(g)=NaN;
                else
                    timeqrs3(g+1:end+1)=timeqrs3(g:end);
                    indexes3(g+1:end+1)= indexes3(g:end);
                    timeqrs3(g)=NaN;
                    indexes3(g)=NaN;
                end
            end
        else
            flag=1;
        end
        if flag==1
            timeqrs=[timeqrs [timeqrs1(g);timeqrs2(g);timeqrs3(g)]]; %#ok<AGROW>
            indexes=[indexes [indexes1(g);indexes2(g);indexes3(g)]]; %#ok<AGROW>
        else
            timeqrs1(g)=NaN;indexes1(g)=NaN;
            timeqrs2(g)=NaN;indexes2(g)=NaN;
            timeqrs3(g)=NaN;indexes3(g)=NaN;
        end
    end
end
lengthaux=([length(timeqrs1) length(timeqrs2) length(timeqrs3)]);
if isempty(g)
    g=0;
end
aux1=find(lengthaux-g >0);
%remains beats in all 3 leads
while length(aux1)==3
%     intervalaux=(g+1):(min(lengthaux(aux1)));
    for g=(g+1):(min(lengthaux(aux1)))       
        [aux auxi]=sort([timeqrs1(g) timeqrs2(g) timeqrs3(g)]);
        if aux(1)>(aux(2)-tol-1)
            if (aux(2)+tol+1)>aux(3)
                timeqrs=[timeqrs [timeqrs1(g);timeqrs2(g);timeqrs3(g)]];  %#ok<AGROW>
                indexes=[indexes [indexes1(g);indexes2(g);indexes3(g)]];  %#ok<AGROW>
            else
                if auxi(3)==1
                    timeqrs1(g+1:end+1)=timeqrs1(g:end);
                    indexes1(g+1:end+1)= indexes1(g:end);
                    indexes1(g)=NaN; 
                    timeqrs1(g)=NaN;
                else if auxi(3)==2
                        timeqrs2(g+1:end+1)=timeqrs2(g:end);
                        indexes2(g+1:end+1)= indexes2(g:end);
                        timeqrs2(g)=NaN;
                        indexes2(g)=NaN;
                    else
                        timeqrs3(g+1:end+1)=timeqrs3(g:end);
                        indexes3(g+1:end+1)= indexes3(g:end);
                        timeqrs3(g)=NaN;
                        indexes3(g)=NaN;
                    end
                end
                timeqrs=[timeqrs [timeqrs1(g);timeqrs2(g);timeqrs3(g)]];  %#ok<AGROW>
                indexes=[indexes [indexes1(g);indexes2(g);indexes3(g)]];  %#ok<AGROW>
            end
        else           
            if auxi(1)==1
                g1=g+1;
                flag=0;
                while  g1<=length(timeqrs1) && flag==0 && timeqrs1(g1)<(aux(2)+tol+1)
                    if (timeqrs1(g1)>(aux(2)-tol-1))
                        flag=1;
                        timeqrs1(g:end)=[timeqrs1(g1:end) NaN*ones(1,(g1-g))];
                        indexes1(g:end)=[indexes1(g1:end) NaN*ones(1,(g1-g))];
                    end
                    g1=g1+1;
                end
                if flag==0       
                    if g1>length(timeqrs1)
                        timeqrs1(g:end)=NaN*ones(1,length(timeqrs1(g:end)));
                        indexes1(g:end)=NaN*ones(1,length(timeqrs1(g:end)));
                    else
                        timeqrs1(g:end)=[timeqrs1((g1-1):end) NaN*ones(1,(g1-1-g))];
                        indexes1(g:end)=[indexes1((g1-1):end) NaN*ones(1,(g1-1-g))];
                        timeqrs1(g)=NaN;
                        indexes1(g)=NaN;
                    end
                end
            elseif auxi(1)==2
                g1=g+1;
                flag=0;
                while  g1<=length(timeqrs2) && flag==0 && timeqrs2(g1)<(aux(2)+tol+1)
                    if (timeqrs2(g1)>(aux(2)-tol-1))
                        flag=1;
                        timeqrs2(g:end)=[timeqrs2(g1:end) NaN*ones(1,(g1-g))];
                        indexes2(g:end)=[indexes2(g1:end) NaN*ones(1,(g1-g))];
                    end
                    g1=g1+1;
                end
                if flag==0       
                    if g1>length(timeqrs1)
                        timeqrs2(g:end)=NaN*ones(1,length(timeqrs2(g:end)));
                        indexes2(g:end)=NaN*ones(1,length(timeqrs2(g:end)));
                    else
                        timeqrs2(g:end)=[timeqrs2((g1-1):end) NaN*ones(1,(g1-1-g))];
                        indexes2(g:end)=[indexes2((g1-1):end) NaN*ones(1,(g1-1-g))];
                        timeqrs2(g)=NaN;
                        indexes2(g)=NaN;
                    end
                end          
            else
                g1=g+1;
                flag=0;
                while  g1<=length(timeqrs3) && flag==0 && timeqrs3(g1)<(aux(2)+tol+1)
                    if (timeqrs3(g1)>(aux(2)-tol-1))
                        flag=1;
                        timeqrs3(g:end)=[timeqrs3(g1:end) NaN*ones(1,(g1-g))];
                        indexes3(g:end)=[indexes3(g1:end) NaN*ones(1,(g1-g))];
                    end
                    g1=g1+1;
                end
                if flag==0       
                    if g1>length(timeqrs1)
                        timeqrs3(g:end)=NaN*ones(1,length(timeqrs3(g:end)));
                        indexes3(g:end)=NaN*ones(1,length(timeqrs3(g:end)));
                    else
                        timeqrs3(g:end)=[timeqrs3((g1-1):end) NaN*ones(1,(g1-1-g))];
                        indexes3(g:end)=[indexes3((g1-1):end) NaN*ones(1,(g1-1-g))];
                        timeqrs3(g)=NaN;
                        indexes3(g)=NaN;
                    end
                end
            end        
            if aux(3)-aux(2)>(tol+1)
                if auxi(3)==1
                    timeqrs1(g+1:end+1)=timeqrs1(g:end);
                    indexes1(g+1:end+1)= indexes1(g:end);
                    indexes1(g)=NaN; 
                    timeqrs1(g)=NaN;
                else if auxi(3)==2
                        timeqrs2(g+1:end+1)=timeqrs2(g:end);
                        indexes2(g+1:end+1)= indexes2(g:end);
                        timeqrs2(g)=NaN;
                        indexes2(g)=NaN;
                    else
                        timeqrs3(g+1:end+1)=timeqrs3(g:end);
                        indexes3(g+1:end+1)= indexes3(g:end);
                        timeqrs3(g)=NaN;
                        indexes3(g)=NaN;
                    end
                end
            else
                flag=1;
            end
            if flag==1
                timeqrs=[timeqrs [timeqrs1(g);timeqrs2(g);timeqrs3(g)]];  %#ok<AGROW>
                indexes=[indexes [indexes1(g);indexes2(g);indexes3(g)]];  %#ok<AGROW>
            else 
                timeqrs1(g)=NaN;indexes1(g)=NaN;
                timeqrs2(g)=NaN;indexes2(g)=NaN;
                timeqrs3(g)=NaN;indexes3(g)=NaN;
            end
        end
    end   
    lengthaux=([length(timeqrs1) length(timeqrs2) length(timeqrs3)]);
    aux1=find(lengthaux-g >0);
end
[min1i,min1]=min(find(~isnan(timeqrs1(end:-1:1)))); %#ok<NASGU,MXFND>
[min2i,min3]=min(find(~isnan(timeqrs2(end:-1:1)))); %#ok<NASGU,MXFND>
[min3i,min3]=min(find(~isnan(timeqrs3(end:-1:1)))); %#ok<NASGU,MXFND>
timeqrs1(length(timeqrs1)-min1i+2:end)=[];
timeqrs2(length(timeqrs2)-min2i+2:end)=[];
timeqrs3(length(timeqrs3)-min3i+2:end)=[];
lengthaux=([length(timeqrs1) length(timeqrs2) length(timeqrs3)]);
aux1=find(lengthaux-g >0);
%remaining beats in 2 leads
if length(aux1)>1
    for intervalauxi=(g+1):(min(lengthaux(aux1)));
        if aux1(1)==1 
            if aux1(2)==2
                [aux auxi]=sort([timeqrs1(intervalauxi) timeqrs2(intervalauxi)]);
                if aux(1)>(aux(2)-tol-1)
                    timeqrs=[timeqrs [timeqrs1(intervalauxi);timeqrs2(intervalauxi);NaN]];  %#ok<AGROW>
                    indexes=[indexes [indexes1(intervalauxi);indexes2(intervalauxi);NaN]];  %#ok<AGROW>
                else
                    if auxi(2)==2
                        timeqrs2(intervalauxi+1:end+1)=timeqrs2(intervalauxi:end);
                        indexes2(intervalauxi+1:end+1)= indexes2(intervalauxi:end);
                        timeqrs2(intervalauxi)=NaN;
                        indexes2(intervalauxi)=NaN;
                    elseif auxi(2)==1
                        timeqrs1(intervalauxi+1:end+1)=timeqrs1(intervalauxi:end);
                        indexes1(intervalauxi+1:end+1)= indexes1(intervalauxi:end);
                        timeqrs1(intervalauxi)=NaN;
                        indexes1(intervalauxi)=NaN;
                    end
                end
            elseif aux1(2)==3
                [aux auxi]=sort([timeqrs1(intervalauxi) timeqrs3(intervalauxi)]);
                if aux(1)>(aux(2)-tol-1)
                    timeqrs=[timeqrs [timeqrs1(intervalauxi);NaN;timeqrs3(intervalauxi)]];  %#ok<AGROW>
                    indexes=[indexes [indexes1(intervalauxi);NaN;indexes3(intervalauxi)]];  %#ok<AGROW>
                else
                    if auxi(2)==2 %25SET08
                        timeqrs3(intervalauxi+1:end+1)=timeqrs3(intervalauxi:end);
                        indexes3(intervalauxi+1:end+1)= indexes3(intervalauxi:end);
                        timeqrs3(intervalauxi)=NaN;
                        indexes3(intervalauxi)=NaN;
                    elseif auxi(2)==1
                        timeqrs1(intervalauxi+1:end+1)=timeqrs1(intervalauxi:end);
                        indexes1(intervalauxi+1:end+1)= indexes1(intervalauxi:end);
                        timeqrs1(intervalauxi)=NaN;
                        indexes1(intervalauxi)=NaN;
                    end
                end
            end
        elseif aux1(1)==2 
            [aux auxi]=sort([timeqrs2(intervalauxi) timeqrs3(intervalauxi)]);
            if aux(1)>(aux(2)-tol-1)
                timeqrs=[timeqrs [NaN;timeqrs2(intervalauxi);timeqrs3(intervalauxi)]];  %#ok<AGROW>
                indexes=[indexes [NaN;indexes2(intervalauxi);indexes3(intervalauxi)]];  %#ok<AGROW>
            else
                if auxi(2)==1%25SET08
                    timeqrs2(intervalauxi+1:end+1)=timeqrs2(intervalauxi:end);
                    indexes2(intervalauxi+1:end+1)= indexes2(intervalauxi:end);
                    timeqrs2(intervalauxi)=NaN;
                    indexes2(intervalauxi)=NaN;
                elseif auxi(2)==2%25SET08
                    timeqrs3(intervalauxi+1:end+1)=timeqrs3(intervalauxi:end);
                    indexes3(intervalauxi+1:end+1)= indexes3(intervalauxi:end);
                    timeqrs3(intervalauxi)=NaN;
                    indexes3(intervalauxi)=NaN;
                end
            end
        end
    end
    lengthaux=([length(timeqrs1) length(timeqrs2) length(timeqrs3)]);
    if isempty(intervalauxi);
        intervalauxi=g;
    end
    aux1=find(lengthaux-intervalauxi >0);
    %remains beats in 2 leads
    while length(aux1)==2
        intervalauxi=intervalauxi+1;
        if aux1(1)==1 
            if aux1(2)==2
                [aux auxi]=sort([timeqrs1(intervalauxi) timeqrs2(intervalauxi)]);
                if aux(1)>(aux(2)-tol-1)
                    timeqrs=[timeqrs [timeqrs1(intervalauxi);timeqrs2(intervalauxi);NaN]];  %#ok<AGROW>
                    indexes=[indexes [indexes1(intervalauxi);indexes2(intervalauxi);NaN]];  %#ok<AGROW>
                else
                    if auxi(2)==2
                        timeqrs2(intervalauxi+1:end+1)=timeqrs2(intervalauxi:end);
                        indexes2(intervalauxi+1:end+1)= indexes2(intervalauxi:end);
                        timeqrs2(intervalauxi)=NaN;
                        indexes2(intervalauxi)=NaN;
                    elseif auxi(2)==1
                        timeqrs1(intervalauxi+1:end+1)=timeqrs1(intervalauxi:end);
                        indexes1(intervalauxi+1:end+1)= indexes1(intervalauxi:end);
                        timeqrs1(intervalauxi)=NaN;
                        indexes1(intervalauxi)=NaN;
                    end
                end
            elseif aux1(2)==3
                [aux auxi]=sort([timeqrs1(intervalauxi) timeqrs3(intervalauxi)]);
                if aux(1)>(aux(2)-tol-1)
                    timeqrs=[timeqrs [timeqrs1(intervalauxi);NaN;timeqrs3(intervalauxi)]];  %#ok<AGROW>
                    indexes=[indexes [indexes1(intervalauxi);NaN;indexes3(intervalauxi)]];  %#ok<AGROW>
                else
                    if auxi(2)==3
                        timeqrs3(intervalauxi+1:end+1)=timeqrs3(intervalauxi:end);
                        indexes3(intervalauxi+1:end+1)= indexes3(intervalauxi:end);
                        timeqrs3(intervalauxi)=NaN;
                        indexes3(intervalauxi)=NaN;
                    elseif auxi(2)==1
                        timeqrs1(intervalauxi+1:end+1)=timeqrs1(intervalauxi:end);
                        indexes1(intervalauxi+1:end+1)= indexes1(intervalauxi:end);
                        timeqrs1(intervalauxi)=NaN;
                        indexes1(intervalauxi)=NaN;
                    end
                end
            end
        elseif aux1(1)==2 
            [aux auxi]=sort([timeqrs2(intervalauxi) timeqrs3(intervalauxi)]);
            if aux(1)>(aux(2)-tol-1)
                timeqrs=[timeqrs [NaN; timeqrs2(intervalauxi);timeqrs3(intervalauxi)]];  %#ok<AGROW>
                indexes=[indexes [NaN; indexes2(intervalauxi);indexes3(intervalauxi)]];  %#ok<AGROW>
            else
                if auxi(2)==2
                    timeqrs2(intervalauxi+1:end+1)=timeqrs2(intervalauxi:end);
                    indexes2(intervalauxi+1:end+1)= indexes2(intervalauxi:end);
                    timeqrs2(intervalauxi)=NaN;
                    indexes2(intervalauxi)=NaN;
                elseif auxi(2)==3
                    timeqrs3(intervalauxi+1:end+1)=timeqrs3(intervalauxi:end);
                    indexes3(intervalauxi+1:end+1)= indexes3(intervalauxi:end);
                    timeqrs3(intervalauxi)=NaN;
                    indexes3(intervalauxi)=NaN;
                end
            end
        end
        lengthaux=([length(timeqrs1) length(timeqrs2) length(timeqrs3)]);
        aux1=find(lengthaux-intervalauxi >0);
    end
    lengthaux=([length(timeqrs1) length(timeqrs2) length(timeqrs3)]);
    g=intervalauxi;
    aux1=find(lengthaux-g >0);   
end
%NOTE THAT ONLY ALIGNED ONES ARE CONSIDERED!!!!!!!!!!!!!!!!!!
intervalaux=(g+1):(min(lengthaux(aux1)));
if length(aux1)>1
    timeqrs1((min(lengthaux(aux1)))+1:end)=[];
    timeqrs2((min(lengthaux(aux1)))+1:end)=[];
    timeqrs3((min(lengthaux(aux1)))+1:end)=[];
    indexes1((min(lengthaux(aux1)))+1:end)=[];
    indexes2((min(lengthaux(aux1)))+1:end)=[];
    indexes3((min(lengthaux(aux1)))+1:end)=[];
    if ~ismember(1,aux1)
        timeqrs1(intervalaux)=NaN*ones(size(intervalaux));
        indexes1(intervalaux)=NaN*ones(size(intervalaux));
        aux=[timeqrs2(intervalaux); timeqrs3(intervalaux) ];
    end    
    if ~ismember(2,aux1)
        timeqrs2(intervalaux)=NaN*ones(size(intervalaux));
        indexes2(intervalaux)=NaN*ones(size(intervalaux));
        aux=[timeqrs1(intervalaux); timeqrs3(intervalaux) ];
    end
    if ~ismember(3,aux1)
        timeqrs3(intervalaux)=NaN*ones(size(intervalaux));
        indexes3(intervalaux)=NaN*ones(size(intervalaux));
        aux=[timeqrs1(intervalaux); timeqrs2(intervalaux) ];
    end  
    timeqrs=[timeqrs [timeqrs1(intervalaux(abs(diff(aux))<(tol+1)));timeqrs2(intervalaux(abs(diff(aux))<(tol+1)));timeqrs3(intervalaux(abs(diff(aux))<(tol+1)))]];
    indexes=[indexes [indexes1(intervalaux(abs(diff(aux))<(tol+1)));indexes2(intervalaux(abs(diff(aux))<(tol+1)));indexes3(intervalaux(abs(diff(aux))<(tol+1)))]];
%     indexes1(intervalaux(abs(diff(aux))>=(tol+1) | isnan(diff(aux))))=[];
%     indexes2(intervalaux(abs(diff(aux))>=(tol+1) | isnan(diff(aux))))=[];
%     indexes3(intervalaux(abs(diff(aux))>=(tol+1) | isnan(diff(aux))))=[];
%     timeqrs1(intervalaux(abs(diff(aux))>=(tol+1) | isnan(diff(aux))))=[];
%     timeqrs2(intervalaux(abs(diff(aux))>=(tol+1) | isnan(diff(aux))))=[];
%     timeqrs3(intervalaux(abs(diff(aux))>=(tol+1) | isnan(diff(aux))))=[];
else
%     timeqrs1((min(lengthaux(aux1))):end)=[];
%     timeqrs2((min(lengthaux(aux1))):end)=[];
%     timeqrs3((min(lengthaux(aux1))):end)=[];
%     indexes1((min(lengthaux(aux1))):end)=[];
%     indexes2((min(lengthaux(aux1))):end)=[];
%     indexes3((min(lengthaux(aux1))):end)=[];
end  
if ~isempty(indexes) %17ABRIL08
indexes=indexes(:,~isnan(indexes(1,:))| ~isnan(indexes(2,:))|~isnan(indexes(3,:)));
end
timeqrs(:,(size(indexes,2)+1):end)=[];
