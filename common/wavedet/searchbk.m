function [new,messages] = searchbk (w, i, thres,freq,thres_alpha,messages)
%performs the searching back when no qrs has been detected for
% 1.5 times the RR before.  i is the beginning sample.  thres
% is the threshold for scale 4;
% Last update: Rute Almeida 27Jan2012

if nargin<6 || (isempty(freq) && (~isfield(messages.setup.wavedet,'freq') || ~(messages.setup.wavedet.freq>0)))
    messages.errors=[messages.errors {'Fatal error in wavedet_3D: no sampling frequency variable found.'}];
    warning(char(messages.errors(end)))
    messages.errors_desc=[messages.errors_desc {'No sampling frequency variable give to seachbk from fiducialf.m.'}];
    new=0;
    return
elseif ~isempty(freq) && isfield(messages.setup.wavedet,'freq') && messages.setup.wavedet.freq>0 && freq~=messages.setup.wavedet.freq
    messages.warnings=[messages.warnings {'sampling frequency used in searchbk is different from messages.setup.wavedet.freq.'}];
elseif isempty(freq)
    freq=messages.setup.wavedet.freq;
end

if ~isfield(messages.setup.wavedet,'refrper') || ~isfield(messages.setup.wavedet,'peakcriteria') || ~isfield(messages.setup.wavedet,'timelapthr') || ~isfield(messages.setup.wavedet,'timelapthr') || ~isfield(messages.setup.wavedet,'intvlthr1') || ~isfield(messages.setup.wavedet,'intvlthr2') || ~isfield(messages.setup.wavedet,'nghbhd')
    messages.errors=[messages.errors {'Fatal error in wavedet_3D: missing parameters on seachbk.'}];
    warning(char(messages.errors(end)))
    messages.errors_desc=[messages.errors_desc {'Setup value is missing on seachbk from fiducialf.m.'}];
    new=0;
    return
end
if ~isfield(messages.setup.wavedet,'intvlthr1_2_sbk')
    messages.setup.wavedet.intvlthr1_2_sbk=2;%reduction on threshold interval for considering redundancy if it subsists
end
if ~isfield(messages.setup.wavedet,'thfraction_sbk')
    if isfield(messages.setup.wavedet,'thfraction')
        messages.setup.wavedet.thfraction_sbk=messages.setup.wavedet.thfraction;
    else
        messages.errors=[messages.errors {'Fatal error in wavedet_3D: missing parameters on seachbk.'}];
        warning(char(messages.errors(end)))
        messages.errors_desc=[messages.errors_desc {'Setup value thfraction is missing on seachbk from fiducialf.m.'}];
        new=0;
        return
    end
end

if ~isfield(messages.setup.wavedet,'nghbhd')
    messages.setup.wavedet.nghbhd=0.025; % neighbourwood fro maximum across scales in sec
end

thfraction_sbk=messages.setup.wavedet.thfraction_sbk;
peakcriteria=messages.setup.wavedet.peakcriteria;
timelap = ceil(messages.setup.wavedet.timelapthr*freq);


thres = thres/2;  % When searching back we take half
thres(4) = thres(4)/2; % so new thres(4) = old thres (4)/4

n = modmax(w(:,4),2,thres(4),0);  % Search maximum moduli at scale 4
neighb = ceil(freq*messages.setup.wavedet.nghbhd);        % Neighbourhood = 25 ms

m = zeros(4,length(n));
signo = sign(w(n,4))';
if ~isempty(n),
    m(4,:) = n';
end
new = []; %#ok<NASGU>


for k = 1:size(m,2),
    window = [max(n(k)-neighb,1), min(n(k)+neighb,size(w,1))];
    n3 = modmax(w(window(1):window(2),3),2,thres(3),signo(k));
    n3 =window(1)-1+n3;
    num = length(n3);
    if num>0,
        if num==1,
            m(3,k)= n3;
        elseif num>1
            if length(find(max(abs(w(n3,3)))./abs(w(n3,3))<peakcriteria))==1,
                [aux1,ind]=max(abs(w(n3,3))); %#ok<ASGLU> % greatest modulus
                m(3,k)= n3(ind);
            else                              % minimum distance
                [aux1,ind]=min(abs(m(4,k)-n3)); %#ok<ASGLU>
                m(3,k)= n3(ind);
            end
        end
    end
end
ind = find(m(3,:)==0);
m(:,ind) = [];
signo(ind) = [];

% Search for maximum moduli in the neighborhood at scale 2
for k = 1:size(m,2),
    window = [max(m(3,k)-neighb,1), min(m(3,k)+neighb,size(w,1))];
    n2 =modmax(w(window(1):window(2),2),2,thres(2),signo(k));
    n2 =window(1)-1+n2;
    num = length(n2);
    if num>0,
        if num==1,
            m(2,k)= n2;
        elseif num>1
            if length(find(max(abs(w(n2,2)))./abs(w(n2,2))<peakcriteria))==1,
                [aux1,ind]=max(abs(w(n2,2))); %#ok<ASGLU> % greatest modulus
                m(2,k)= n2(ind);
            else                            % shortest distance
                [aux1,ind]=min(abs(m(3,k)-n2)); %#ok<ASGLU>
                m(2,k)= n2(ind);
            end
        end
    end
end
ind = find(m(2,:)==0);
m(:,ind) = [];
signo(ind) = [];

for k = 1:size(m,2),
    window = [max(m(2,k)-neighb,1), min(m(2,k)+neighb,size(w,1))];
    n1 =modmax(w(window(1):window(2),1),2,thres(1),signo(k));
    n1 =window(1)-1+n1;
    num = length(n1);
    if num>0,
        if num==1,
            m(1,k)= n1;
        elseif num>1
            if length(find(max(abs(w(n1,1)))./abs(w(n1,1))<peakcriteria))==1,
                [aux1,ind]=max(abs(w(n1,1)));%#ok<ASGLU> % greatest modulus
                m(1,k)= n1(ind);
            else                             % shortest distance
                [aux1,ind]=min(abs(m(2,k)-n1)); %#ok<ASGLU>
                m(1,k)= n1(ind);
            end
        end
    end
end

ind = find(m(1,:)==0);               %Discard all maximum lines with no
m(:,ind) = [];                       % associated maximum at scale 1
signo(ind) = [];

% Regularity Exponent Validation
% alpha proportional to log(a3(nk3))-log(a1(nk1))
alpha = log(abs(w(m(3,:),3)));  %%!!!!
% alpha = log(abs(w(m(3,:),3))) - log (abs(w(m(1,:),1)));  !!!
ind = find(alpha <= thres_alpha-thfraction_sbk);   %%%% !!!
m(:,ind) = [];
signo(ind) = [];

thresinterval = ceil(messages.setup.wavedet.intvlthr2 * freq);         % 120 ms. (Li)

if size(m,2)>2,
    ind = find( ((m(1,2:end-1)-m(1,1:end-2))>thresinterval) ...
        &      (m(1,3:end)- m(1,2:end-1))>thresinterval)+1;
    if (m(1,2)-m(1,1))>thresinterval,
        ind = [1 ind];
    end
    if (m(1,end)-m(1,end-1))>thresinterval,
        ind = [ind size(m,2)];
    end
    m(:,ind) =[];                     % Discard isolated maximum lines
    signo(ind) = [];
    
elseif size(m,2)==2,                 % If only two lines
    if (m(1,2)-m(1,1))>thresinterval, % discard them if too separated
        m(:,1:2)=[];
        signo(1:2)=[];
    end
elseif size(m,2)==1,                % If only one, discard it
    m(:,1)=[];
    signo(1)=[];
end

% Threshold interval for considering redundancy
redundant= [];
thresinterval = ceil(messages.setup.wavedet.intvlthr1*freq);           % 120 ms. (li)

for l = find(signo>0),                      % For each positive maximum line
    if ~any(redundant ==l),                   % If it has not been declared redundant yet
        ind=find((m(3,:)>m(3,l)-messages.setup.wavedet.intvlthr1_2_sbk*thresinterval)&(m(3,:)<m(3,l)+messages.setup.wavedet.intvlthr1_2_sbk*thresinterval)&signo>0);
        % index of positive lines near the present one (including it)
        if length(ind)>1,                        % If more than one --> redundancy
            [mx,ind2]= max(abs(w(m(3,ind),3))); %#ok<ASGLU>
            ind(ind2)=[];
            redundant = [redundant ind];             %#ok<AGROW> % All but the greatest are redundant
        end
    end
end

m(:,redundant)=[];                         % Discard redundant lines
signo(redundant)=[];
redundant= [];


for l = find(signo>0),              % For each remaining positive maximum line
    ind = find((m(3,:)>m(3,l)-thresinterval)&(m(3,:)<m(3,l)+thresinterval)&signo<0);
    % Search for negative minima near it
    if length(ind)>1,               % If more than one ---> redundancy
        aux = abs(w(m(3,ind),3)./(m(3,ind)'-m(3,l)));
        % auxiliary variable: height over distance
        [mx,indmx]=max(aux);
        ind2=ind;
        aux(indmx)=[]; ind2(indmx)=[];
        if all((mx./aux)>peakcriteria),        % RULE 2 (see PFC or Li's paper)
            redundant = [redundant ind2]; %#ok<AGROW>
        else
            [aux,aux2]=min(abs(m(3,l)-m(3,ind))); %#ok<ASGLU>
            ind(aux2)=[];              % RULE 1 (see PFC or Li's paper)
            redundant = [redundant ind]; %#ok<AGROW>
        end
    end
end
m(:,redundant)=[];                  % Discard redundant lines
signo(redundant)=[];
redundant = [];
for l = find(signo<0),              % For each remaining negative minimum line
    ind = find((m(3,:)>m(3,l)-thresinterval)&(m(3,:)<m(3,l)+thresinterval)&signo>0);
    % Search for positive maxima near it
    if length(ind)>1,                % If more than one ----> redundancy
        aux = abs(w(m(3,ind),3)./(m(3,ind)'-m(3,l)));
        % auxiliary variable: height over distance
        [mx,indmx]=max(aux);
        ind2=ind;
        aux(indmx)=[]; ind2(indmx)=[];
        if all((mx./aux)>peakcriteria),
            redundant = [redundant ind2]; %#ok<AGROW>   % RULE 2
        else
            [aux,aux2]=min(abs(m(3,l)-m(3,ind))); %#ok<ASGLU>
            ind(aux2)=[];
            redundant = [redundant ind]; %#ok<AGROW>    % RULE 1
        end
    end
end
m(:,redundant)=[];
signo(redundant)=[];

%%%%%%%%%%%%%% isolated maximum lines resulting from Discarded
%%%%%%%%%%%%%% redundant OUT1011

ind = find( ((m(1,2:end-1)-m(1,1:end-2))>thresinterval) ...
    &      (m(1,3:end)- m(1,2:end-1))>thresinterval)+1;


if size(m,2)<2 && ~isempty(ind), ind=1; end  % when only one maximum %16DEZ08
m(:,ind) = [];                     %Discard isolated maximum lines
signo(ind)=[];
if size(m,2)<2
    m=[];
    signo=[];
end

eliminar=[];
%%%%extra protection% OUT2011
for ii=1:size(m,2)
    pa = picant(w(max(1,m(2,ii)-round(messages.setup.wavedet.pictime*freq)):m(2,ii),2),m(2,ii));
    % first peak before detected qrs position at scale 2
    pp = picpost(w(m(2,ii):min(size(w,1),m(2,ii)+round(messages.setup.wavedet.pictime*freq)),2),m(2,ii));
    
    if isempty(pa) || isempty(pp)
        eliminar=[eliminar ii]; %#ok<AGROW>
    end
end

m(:,eliminar)=[];
signo(eliminar)=[];



% QRS peak detection / wavelet Zero cross detection
time = [];
aux=[]; %Rute26Jun09
if length(signo)>1,
    for l = find(signo>0),
        if (l==1)
            if (signo(2)<0)&&(m(1,2)-m(1,1)<thresinterval),
                ind = zerocros(w(m(1,l):min(m(1,l)+timelap,size(w,1)),1));
                % Zero crossing at scale 1
                time = [time ind+m(1,l)-1]; %#ok<AGROW>
                aux=[aux abs(w(m(1,l))-w(m(1,l+1)))];  %#ok<AGROW> %Rute26Jun09
            end
        elseif (l==size(m,2))                % Special case: last line
            if (signo(l-1)<0)&&(m(1,end)-m(1,end-1)<thresinterval),
                ind = zerocros(w(m(1,l-1):min(m(1,l-1)+timelap,size(w,1)),1));
                time = [time ind+m(1,l-1)-1]; %#ok<AGROW>
                aux=[aux abs(w(m(1,l-1))-w(m(1,l)))];  %#ok<AGROW>  %Rute26Jun09
            end
        elseif signo(l+1)<0 && ((signo(l-1)>0) ...
                || ((m(1,l+1)-m(1,l))<(m(1,l)-m(1,l-1))))
            ind = zerocros(w(m(1,l):min(m(1,l)+timelap,size(w,1)),1));
            time = [time ind+m(1,l)-1]; %#ok<AGROW>
            aux=[aux abs(w(m(1,l))-w(m(1,l+1)))]; %#ok<AGROW>  %Rute26Jun09
        elseif signo(l-1)<0,
            ind = zerocros(w(m(1,l-1):min(m(1,l-1)+timelap,size(w,1)),1));
            time = [time ind+m(1,l-1)-1]; %#ok<AGROW>
            aux=[aux abs(w(m(1,l-1))-w(m(1,l)))]; %#ok<AGROW> %Rute 26Jun09
        end
    end
end

% Refractary period after a QRS detection (200 ms).
rr = (time(2:end)-time(1:end-1));
ind = find(rr<ceil(messages.setup.wavedet.refrper(end)*freq));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RUTE 26Jun09
%time(ind+1)=[]; % always the second one is eliminated!!! to be changed
for auxi=1:length(ind) %RUTE 27Jun11
    [M,ii]=min([aux(ind(auxi)) aux(ind(auxi)+1)]); %#ok<ASGLU> % 18JUL2011
    if ii==2, ind(auxi)=ind(auxi)+1; end
end
time(ind)=[]; %RUTE 27Jun11
aux(ind)=[]; %#ok<NASGU>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new = time +i -1;