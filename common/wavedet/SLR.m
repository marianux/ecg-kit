function [marksSLR,ind,messages]=SLR(marks,fa,QRSLIMIT,k,delta,messages)
% post processing rule for delineation
%
% [marksSLR,ind,messages]=SLR(marks,fa,QRSLIMIT,k,delta,messages)
%
% INPUT:
% marks: structure with a matrix of annotations in samples; in each field
%        each matrix of annotations has one line per beat and 10 colunms
%        corresponding respectively to P onset, P peak, P end, QRS onset, QRS
%        main peak, QRS end, T onset, T peak, T prima, T end
% fa: sampling frequency of the recording.
% QRSLIMIT: admited diference between marks of the same beat in samples (optional)
%           by default QRSLIMIT= 250 ms.
% k: number of neighbours to consider in SLR rule (optional)
%    by default k=3 
% delta: neighbourhood to consider in SLR rule in samples (optional)
%    by default delta = [4 nan 4 12 nan 10 12 nan nan 12] ms for 10 marks
%    or  delta = [4 nan nan 4 nan 12  nan nan nan nan 10 12 nan nan 12] ms for 15 marks
% messages:- errors and warnings (optional)
%
% OUTPUT:
% marksSLR: multilead based on marks obtained by SLR rule
% ind: indexes of the initial column vectors for each set, corresponding to the sincronized beats
%      each line of ind corresponds to a beat and each column to a set
%      if a beat is absent of a set
% messages.errors: errors during the processing
% messages.errors_desc:errors description
% messages.status: sucess of the procesing (1) ou fail (0)
% messages.status: parameters considered
% Last update: DEZ2011

marksSLR=[];
ind=[];
if nargin<2
    messages.errors= 'Fatal error in SLR.';
    messages.errors_desc= 'Inputs marks and fa required';
    messages.status=0;
    return
else
    if size(marks{1},2)~=10 && size(marks{1},2)~=15
        messages.errors= 'Fatal error in SLR.';
        messages.errors_desc= 'Input marks should include 10 or 15 marks.';
        messages.status=0;
        return
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
    if ~isfield(messages,'setup')
        messages.setup=[];
    end
    if ~isfield(messages.setup,'SLR')
        messages.setup.SLR=[];
    end
end

if nargin<5
    if size(marks{1},2)==10
        delta = [4 nan 4 12 nan 10 12 nan nan 12]*fa/1000;
    elseif size(marks{1},2)==15
        delta = [4 nan nan 4 nan 12  nan nan nan nan 10 12 nan nan 12]*fa/1000;
    end
    if nargin<4
        k=min(3,size(marks{1},1)-1);
        k=2;
        if nargin<3
            QRSLIMIT=250*fa/1000;
        end
    end
    messages.setup.SLR.QRSlimit=QRSLIMIT;
    messages.setup.SLR.delta=delta;
    messages.setup.SLR.k=k;
    if nargin<3
        if isfield(messages,'setup') & isfield(messages.setup,'SLR')
            if isfield(messages.setup.SLR,'QRSlimit')
                QRSLIMIT=messages.setup.SLR.QRSlimit;
            else
                QRSLIMIT=[];
            end
            if isfield(messages.setup.SLR,'k')
                k=messages.setup.SLR.k;
            else
                k=[];
            end
            if isfield(messages.setup.SLR,'delta')
                if length(delta)==size(marks{1},2)
                    delta= messages.setup.SLR.delta;
                else
                    delta=[];
                end
            end
        else
            QRSLIMIT=[];
            k=[];
            delta=[];
        end
    end
end

if isempty(QRSLIMIT) | ~isnumeric(QRSLIMIT) | isnan(QRSLIMIT)
    QRSLIMIT=250*fa/1000;
    messages.setup.SLR.QRSlimit=QRSLIMIT;
end
if isempty(k) | ~isnumeric(k) | isnan(k)
    k=min(3,size(marks{1},1)-1);
    messages.setup.SLR.k=k;
end
if isempty(delta) | size(delta)~=size(marks{1},2)
    if size(marks{1},2)==10
        delta = [4 nan 4 12 nan 10 12 nan nan 12]*fa/1000;
    elseif size(marks{1},2)==15
        delta = [4 nan nan 4 nan 12  nan nan nan nan 10 12 nan nan 12]*fa/1000;
    end
    messages.setup.SLR.delta=delta;
end
    
if ~isfield(messages.setup.SLR,'refrper')
    if isfield(messages.setup,'wavedet') & isfield(messages.setup.wavedet,'refrper')
        messages.setup.SLR.refrper=messages.setup.wavedet.refrper;
    else
        messages.setup.SLR.refrper=0.275;
    end
end


if ~isfield(messages.setup.SLR,'tolerance')
    messages.setup.SLR.tolerance=10*fa/1000; % 10 msec
end

if isstruct(marks)
    leads = fieldnames(marks);
    messages.setup.nleads=length(leads);
else
   messages.setup.nleads = size(marks,2);
end

if length(delta)==10
    messages.setup.marks_desc={'(','P',')','(','QRS',')','(','T','T''',')'};
%     nmarks=10;
    peaks=[2 5 8 9];
    onsets=[1 4 6];
    ends=[3 6 10];
    refmark=5; %JB 27JUL2011
elseif length(delta)==15
    messages.setup.marks_desc={'(','P','P''',')','QRS_first_peak','(','Q','R','S','R''',')','(','T','T''',')'};
%     nmarks=15;
    peaks=[2 3 8 13 14];
    onsets=[1 6 12];
    ends=[4 11 15];
   % refmark=8; %JB 27JUL2011
   refmark=[8 9 7 10];%RUTE %29JUL2011 porque la onda R puede ni existir en una derivacion especifica!
end
nmarks=length(delta);
beats=0;
% refmark=5;  %JB 27JUL2011
 refmark_n=refmark; %DEZ2011
for g=1:messages.setup.nleads
    if isstruct(marks)
        f = getfield(marks, char(leads(g))); %#ok<GFLD>
    else
        f = marks{g};
    end
    if length(f)> size(beats,1)
        beats((size(beats,1)+1):length(f),1:g)=zeros(length((size(beats,1)+1):length(f)),g);
    end
    % beats(1:size(f,1),g)=f(:,refmark(1));  %19ABRIL2010
    for n=1:size(f,1)%DEZ2011
        refmark(n)=refmark_n(find(~isnan(f(n,refmark_n)),1,'first'));%DEZ2011
        beats(n,g)=f(n,refmark(n)); %29JUL2011 %DEZ2011
    end
  %beats(1:size(f,1),g)=f(:,refmark); %29JUL2011     
end
beats(beats==0)=NaN;
[ind,messages1]=sincronizabeats(beats,QRSLIMIT); %27JUL2011
ind(sum(~isnan(ind),2)<k+1,:) = [];
% nbeats=size(ind,1);

%JB 27JUL2011
beats=NaN*ones(size(ind));
for g=1:size(ind,2)
    if isstruct(marks)
        f = getfield(marks, char(leads(g))); %#ok<GFLD>
    else
        f = marks{g};
    end
    beats(~isnan(ind(:,g)),g)=f(ind(~isnan(ind(:,g)),g),5);
end
%ind(diff(nanmedian(beats,2)) < 200,:) = [];%JB 27JUL2011

ind(find(diff(nanmedian(beats,2)) < messages.setup.SLR.refrper*fa)+1,:) = [];%RUTE 10AGO 2011


nbeats=size(ind,1);
marksSLR=NaN.*ones(nbeats,nmarks);
for mark=peaks
    beats=NaN*ones(size(ind));
    for g=1:size(ind,2)
        if isstruct(marks)
            f = getfield(marks, char(leads(g))); %#ok<GFLD>
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
            f = getfield(marks, char(leads(g))); %#ok<GFLD>
        else
            f = marks{g};
        end
        QRS{1}(~isnan(ind(:,g)),g)=f(ind(~isnan(ind(:,g)),g),mark(1)); 
        QRS{2}(~isnan(ind(:,g)),g)=f(ind(~isnan(ind(:,g)),g),mark(2)); 
        QRS{3}(~isnan(ind(:,g)),g)=f(ind(~isnan(ind(:,g)),g),mark(3)); 
        QRS{4}(~isnan(ind(:,g)),g)=f(ind(~isnan(ind(:,g)),g),mark(4)); 
    end
    

    for beat=1:nbeats
        AA=[QRS{1}(beat,:);QRS{2}(beat,:);QRS{3}(beat,:);QRS{4}(beat,:)];
        %look for small R waves that should go with Q waves
        %RS that is an inverted QR or RSR' that is an inverted QRS
        ii=find(AA(2,:)-nanmedian(AA(2,:))<-messages.setup.SLR.tolerance);
        QRS{1}(beat,ii)=QRS{2}(beat,ii); 
        QRS{2}(beat,ii)=QRS{3}(beat,ii); 
        QRS{3}(beat,ii)=QRS{4}(beat,ii); 
        QRS{4}(beat,ii)=NaN; 
        AA=[QRS{1}(beat,:);QRS{2}(beat,:);QRS{3}(beat,:);QRS{4}(beat,:)];

        %look for Q waves that should go with R waves
        %QS complexes or QR that is an inverted RS
        ii=find(AA(1,:)-nanmedian(AA(1,:))>messages.setup.SLR.tolerance);
        QRS{4}(beat,ii)=QRS{3}(beat,ii); 
        QRS{3}(beat,ii)=QRS{2}(beat,ii); 
        QRS{2}(beat,ii)=QRS{1}(beat,ii); 
        QRS{1}(beat,ii)=NaN; 
        AA=[QRS{1}(beat,:);QRS{2}(beat,:);QRS{3}(beat,:);QRS{4}(beat,:)];

        %look for R' waves that should go with S waves
        ii=find(AA(4,:)-nanmedian(AA(4,:))<-messages.setup.SLR.tolerance);
        QRS{1}(beat,ii)=QRS{2}(beat,ii); 
        QRS{2}(beat,ii)=QRS{3}(beat,ii); 
        QRS{3}(beat,ii)=QRS{4}(beat,ii); 
        QRS{4}(beat,ii)=NaN;  
        AA=[QRS{1}(beat,:);QRS{2}(beat,:);QRS{3}(beat,:);QRS{4}(beat,:)];
        marksSLR(beat,mark)=nanmedian(AA');
    end
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
            kk =sort(kk); % sorted the mark found in the leads in beat
            while (length(kk)>k && kk(k+1)>(kk(1)+delta(mark))), %eliminate the first mark until kk satisfy the rule or length(kk)<k
                kk(1)=[];
            end
            %if length(kk)<k there is not a multilead mark
            if length(kk)==k
                marksSLR(beat,mark)=nan;
            else
                marksSLR(beat,mark)=kk(1);
            end
        else
            marksSLR(beat,mark)=nan;
        end
    end

end
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
            kk =flipud(sort(kk)')'; %Rute 19AB2010
            while (length(kk)>k && kk(k+1)<(kk(1)-delta(mark))),
                kk(1)=[];
            end
            if length(kk)==k
                marksSLR(beat,mark)=nan;
            else
                marksSLR(beat,mark)=kk(1);
            end
        else
            marksSLR(beat,mark)=nan;
        end
    end
end