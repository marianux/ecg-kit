function [positionSLR,marksSLR,ind,messages]=SLR2(marks,messages)
% post processing rule for delineation
%
% [positionSLR,marksSLR,ind,messages]=SLR2(marks,messages)
%
% INPUT:
% marks: structure with a matrix of annotations in samples; in each field
%        each matrix of annotations has one line per beat and 10 colunms
%        corresponding respectively to P onset, P peak, P end, QRS onset, QRS
%        main peak, QRS end, T onset, T peak, T prima, T end
% messages.setup.SLR.freq: sampling frequency of the recording.
% messages.setup.SLR.QRSlimit: admited diference between marks of the same beat in ms (optional)
%           by default QRSlimit= 250 ms.
% messages.setup.SLR.k: number of neighbours to consider in SLR rule (optional)
%    by default k=3 
% messages.setup.SLR.delta: neighbourhood to consider in SLR rule in ms (optional)
%    by default delta = [4 nan 4 12 nan 10 12 nan nan 12] ms for 10 marks
%    or  delta = [4 nan nan 4 nan 12  nan nan nan nan 10 12 nan nan 12] ms for 15 marks
% messages.setup.wavedet.refrper: minimum time between two admisible beats, by default 257ms
% messages.setup.SLR.tolerance = time in ms to define each individual QRS wave, by default 10ms
% 
% OUTPUT:
%positionSLR: struture with multilead based on marks obtained by SLR rule
% marksSLR: matrix multilead based on marks obtained by SLR rule
% ind: indexes of the initial column vectors for each set, corresponding to the sincronized beats
%      each line of ind corresponds to a beat and each column to a set
%      if a beat is absent of a set
% messages.errors: errors during the processing
% messages.errors_desc:errors description
% messages.status: sucess of the procesing (1) ou fail (0)
% messages.status: parameters considered
% Last update: 27JAN2012
positionSLR=[];
marksSLR=[];
ind=[];
if nargin<2
    messages.errors= 'Fatal error in SLR.';
    messages.errors_desc= 'Inputs marks and messages.setup.SLR.freq required';
    messages.status=0;
    return
else
    if size(marks{1},2)~=10 && size(marks{1},2)~=15
        messages.errors= 'Fatal error in SLR.';
        messages.errors_desc= 'Input marks should include 10 or 15 marks.';
        messages.status=0;
        return
    end
    if ~isfield(messages,'setup') || ((~isfield(messages.setup,'SLR') ||...
            (~isfield(messages.setup.SLR,'freq') && ~isfield(messages.setup.SLR,'fa')))  &&...
            (~isfield(messages.setup.wavedet,'wavedet') || ( isfield(messages.setup.wavedet,'wavedet') &&...
            ~isfield(messages.setup.wavedet,'freq')))) %
        messages.errors= 'Fatal error in SLR.';
        messages.errors_desc= 'Inputs messages.setup.SLR.freq required';
        messages.status=0;
        return
    else
        try
            fa=messages.setup.SLR.freq;
        catch me
            me.message = 'Not exist messages.setup.SLR.freq';
            try
                fa=messages.setup.SLR.fa;
            catch me
                me.message = 'Not exist messages.setup.SLR.fa';
                fa= messages.setup.wavedet.freq;
            end
        end
    end
end
messages.status=1;
if ~isfield(messages,'errors')
    messages.errors=[];
end
if ~isfield(messages,'errors_desc')
    messages.errors_desc=[];
end
if ~isfield(messages,'warnings')
    messages.warnings=[];
end
if ~isfield(messages,'status')
    messages.status=1;
end
if ~isfield(messages.setup.SLR,'k')
    messages.setup.SLR.k=min(3,size(marks{1},1)-1);
end
if ~isfield(messages.setup.SLR,'QRSlimit')
    messages.setup.SLR.QRSlimit=250;
end
if ~isfield(messages.setup.SLR,'delta')
    if size(marks{1},2)==10
        messages.setup.SLR.delta = [4 nan 4 12 nan 10 12 nan nan 12];
    elseif size(marks{1},2)==15
        messages.setup.SLR.delta = [4 nan nan 4 nan 12  nan nan nan nan 10 12 nan nan 12];
    end
end

QRSlimit= messages.setup.SLR.QRSlimit*fa/1000;
delta= messages.setup.SLR.delta*fa/1000;
k= messages.setup.SLR.k;
if isempty(k) || ~isnumeric(k) || isnan(k)
    k=min(3,size(marks{1},1)-1);
    messages.setup.SLR.k=k;
end
if isempty(QRSlimit) || ~isnumeric(QRSlimit) || isnan(QRSlimit)
    messages.setup.SLR.QRSlimit=250;
    QRSlimit=messages.setup.SLR.QRSlimit*fa/1000;
end
if isempty(delta) | size(delta)~=size(marks{1},2) %#ok<OR2>
    if size(marks{1},2)==10
        messages.setup.SLR.delta = [4 nan 4 12 nan 10 12 nan nan 12];
    elseif size(marks{1},2)==15
        messages.setup.SLR.delta = [4 nan nan 4 nan 12  nan nan nan nan 10 12 nan nan 12];
    end
    delta= messages.setup.SLR.delta*fa/1000;
end
    
if ~isfield(messages.setup.SLR,'refrper')
    if isfield(messages.setup,'wavedet') && isfield(messages.setup.wavedet,'refrper')
        messages.setup.SLR.refrper=messages.setup.wavedet.refrper*fa/1000;
    else
        messages.setup.SLR.refrper=275*fa/1000;
    end
end
if ~isfield(messages.setup.SLR,'tolerance')
    messages.setup.SLR.tolerance=10*fa/1000;
end

if isstruct(marks)
    leads = fieldnames(marks);
    messages.setup.SLR.nleads=length(leads);
else
   messages.setup.SLR.nleads = size(marks,2);
end

if length(delta)==10
    messages.setup.SLR.marks_desc={'(','P',')','(','QRS',')','(','T','T''',')'};
    peaks=[2 5 8 9];
    onsets=[1 4 6];
    ends=[3 6 10];
    refmark=5; 
elseif length(delta)==15
    messages.setup.SLR.marks_desc={'(','P','P''',')','QRS_first_peak','(','Q','R','S','R''',')','(','T','T''',')'};
    peaks=[2 3 8 13 14];
    onsets=[1 6 12];
    ends=[4 11 15];
    refmark=[8 9 7 10];%RUTE % porque la onda R puede ni existir en una derivacion especifica!
end
nmarks=length(delta);
beats=0;
refmark_n=refmark;
 
for g=1:messages.setup.SLR.nleads
    if isstruct(marks)
        f = getfield(marks, char(leads(g)));  %#ok<GFLD>
    else
        f = marks{g};
    end
    if length(f)> size(beats,1)
        beats((size(beats,1)+1):length(f),1:g)=zeros(length((size(beats,1)+1):length(f)),g);
    end
    for n=1:size(f,1)
        refmark(n)=refmark_n(find(~isnan(f(n,refmark_n)),1,'first')); %#ok<AGROW>
        beats(n,g)=f(n,refmark(n));  %#ok<AGROW>
    end     
end

beats(beats==0)=NaN;
[ind,messages1]=sincronizabeats(beats,QRSlimit); 
if messages1.status == 0
    return;
end
ind(sum(~isnan(ind),2)<k+1,:) = [];

beats=NaN*ones(size(ind));
for g=1:size(ind,2)
    if isstruct(marks)
        f = getfield(marks, char(leads(g)));  %#ok<GFLD>
    else
        f = marks{g};
    end
    beats(~isnan(ind(:,g)),g)=f(ind(~isnan(ind(:,g)),g),5);
end
ind(find(diff(nanmedian(beats,2)) < messages.setup.SLR.refrper)+1,:) = [];

nbeats=size(ind,1);
marksSLR=NaN.*ones(nbeats,nmarks);
marksSLR_lead=NaN.*ones(nbeats,nmarks);
for mark=peaks
    beats=NaN*ones(size(ind));
    for g=1:size(ind,2)
        if isstruct(marks)
            f = getfield(marks, char(leads(g)));  %#ok<GFLD>
        else
            f = marks{g};
        end
        beats(~isnan(ind(:,g)),g)=f(ind(~isnan(ind(:,g)),g),mark);
    end
    marksSLR(:,mark)=nanmedian(beats,2);
end

if  nmarks==15;
    mark=7:10;
    QRS{1}=NaN*ones(size(ind));
    QRS{2}=QRS{1};
    QRS{3}=QRS{1};
    QRS{4}=QRS{1};
    for g=1:size(ind,2)
        if isstruct(marks)
            f = getfield(marks, char(leads(g)));  %#ok<GFLD>
        else
            f = marks{g};
        end
        QRS{1}(~isnan(ind(:,g)),g)=f(ind(~isnan(ind(:,g)),g),mark(1));  %#ok<AGROW>
        QRS{2}(~isnan(ind(:,g)),g)=f(ind(~isnan(ind(:,g)),g),mark(2)); %#ok<AGROW>
        QRS{3}(~isnan(ind(:,g)),g)=f(ind(~isnan(ind(:,g)),g),mark(3)); %#ok<AGROW>
        QRS{4}(~isnan(ind(:,g)),g)=f(ind(~isnan(ind(:,g)),g),mark(4)); %#ok<AGROW>
    end  
    for beat=1:nbeats
        AA=[QRS{1}(beat,:);QRS{2}(beat,:);QRS{3}(beat,:);QRS{4}(beat,:)];
        %look for small R waves that should go with Q waves
        %RS that is an inverted QR or RSR' that is an inverted QRS
        ii=find(AA(2,:)-nanmedian(AA(2,:))<-messages.setup.SLR.tolerance);
        QRS{1}(beat,ii)=QRS{2}(beat,ii); %#ok<AGROW>
        QRS{2}(beat,ii)=QRS{3}(beat,ii); %#ok<AGROW>
        QRS{3}(beat,ii)=QRS{4}(beat,ii); %#ok<AGROW>
        QRS{4}(beat,ii)=NaN; %#ok<AGROW>
        AA=[QRS{1}(beat,:);QRS{2}(beat,:);QRS{3}(beat,:);QRS{4}(beat,:)];
        %look for Q waves that should go with R waves
        %QS complexes or QR that is an inverted RS
        ii=find(AA(1,:)-nanmedian(AA(1,:))>messages.setup.SLR.tolerance);
        QRS{4}(beat,ii)=QRS{3}(beat,ii); %#ok<AGROW>
        QRS{3}(beat,ii)=QRS{2}(beat,ii); %#ok<AGROW>
        QRS{2}(beat,ii)=QRS{1}(beat,ii); %#ok<AGROW>
        QRS{1}(beat,ii)=NaN; %#ok<AGROW>
        AA=[QRS{1}(beat,:);QRS{2}(beat,:);QRS{3}(beat,:);QRS{4}(beat,:)];
        %look for R' waves that should go with S waves
        ii=find(AA(4,:)-nanmedian(AA(4,:))<-messages.setup.SLR.tolerance);
        QRS{1}(beat,ii)=QRS{2}(beat,ii); %#ok<AGROW>
        QRS{2}(beat,ii)=QRS{3}(beat,ii); %#ok<AGROW>
        QRS{3}(beat,ii)=QRS{4}(beat,ii); %#ok<AGROW>
        QRS{4}(beat,ii)=NaN;  %#ok<AGROW>
        AA=[QRS{1}(beat,:);QRS{2}(beat,:);QRS{3}(beat,:);QRS{4}(beat,:)];
        marksSLR(beat,mark)=nanmedian(AA,2)';
    end
end
if size(marksSLR,2) == 15
    positionSLR.P = marksSLR(:,2);
    positionSLR.Pprima = marksSLR(:,3);
    positionSLR.qrs = marksSLR(:,8);
    positionSLR.QRSon = marksSLR(:,6);
    positionSLR.Q = marksSLR(:,7);
    positionSLR.R = marksSLR(:,8);
    positionSLR.S = marksSLR(:,9);
    positionSLR.Rprima = marksSLR(:,10);
    positionSLR.T = marksSLR(:,13);
    positionSLR.Tprima = marksSLR(:,14);
else
    positionSLR.P = marksSLR(:,2);
    positionSLR.qrs = marksSLR(:,5);
    positionSLR.T = marksSLR(:,8);
    positionSLR.Tprima = marksSLR(:,9);
end
for mark=onsets, %onsets
    beats=NaN.*ones(nbeats,nmarks);
    for g=1:size(ind,2)
        if isstruct(marks)
            f = getfield(marks, char(leads(g))); %#ok<GFLD>
        else
            f = marks{g};
        end
        beats(~isnan(ind(:,g)),g)=f(ind(~isnan(ind(:,g)),g),mark);
    end
    aux= sum(~isnan(beats),2);
    for beat=1:nbeats
        if aux(beat)>k, % if the mark was found in more than k of the leads
            kk = beats(beat,:);
            kk(isnan(kk))=[];
            [kk kki] =sort(kk); % sorted the mark found in the leads in beat
            while (length(kk)>k && kk(k+1)>(kk(1)+delta(mark))), %eliminate the first mark until kk satisfy the rule or length(kk)<k
                kk(1)=[];
                kki(1)=[];
            end
            %if length(kk)<k there is not a multilead mark
            if length(kk)==k
                marksSLR(beat,mark)=nan;
            else
                marksSLR(beat,mark)=kk(1);
                marksSLR_lead(beat,mark)=kki(1);
            end
        else
            marksSLR(beat,mark)=nan;
        end
    end
end
positionSLR.Pon = marksSLR(:,onsets(1));
positionSLR.QRSon = marksSLR(:,onsets(2));
positionSLR.Ton = marksSLR(:,onsets(3));
positionSLR.Pon_lead = marksSLR_lead(:,onsets(1));
positionSLR.QRSon_lead = marksSLR_lead(:,onsets(2));
positionSLR.Ton_lead = marksSLR_lead(:,onsets(3));

for mark=ends, % ends
    beats=NaN.*ones(nbeats,nmarks);
    for g=1:size(ind,2)
        if isstruct(marks)
            f = getfield(marks, char(leads(g))); %#ok<GFLD>
        else
            f = marks{g};
        end
        beats(~isnan(ind(:,g)),g)=f(ind(~isnan(ind(:,g)),g),mark);
    end
    aux= sum(~isnan(beats),2);
    for beat=1:nbeats
        if aux(beat)>k, %
            kk = beats(beat,:);
            kk(isnan(kk))=[];
            [kk kki] =sort(kk,'descend');
           % kk =flipud(sort(kk)')'; 
            while (length(kk)>k && kk(k+1)<(kk(1)-delta(mark))),
                kk(1)=[];
                kki(1)=[];
            end
            if length(kk)==k
                marksSLR(beat,mark)=nan;
            else
                marksSLR(beat,mark)=kk(1);
                marksSLR_lead(beat,mark)=kki(1);
            end
        else
            marksSLR(beat,mark)=nan;
        end
    end
end
positionSLR.Poff = marksSLR(:,ends(1));
positionSLR.QRSoff = marksSLR(:,ends(2));
positionSLR.Toff = marksSLR(:,ends(3));
positionSLR.Poff_lead = marksSLR_lead(:,ends(1));
positionSLR.QRSoff_lead = marksSLR_lead(:,ends(2));
positionSLR.Toff_lead = marksSLR_lead(:,ends(3));
