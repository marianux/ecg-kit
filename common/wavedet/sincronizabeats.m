function [ind,messages]=sincronizabeats(beats, QRSLIMIT,messages)
% This function sincronizes a set of marks from diferent delineation
% results, assuming that marks difereng less than QRSLIMIT samples correspond to
% the same beat
%
% ind=sincronizabeats(beats, QRSLIMIT)
%
% INPUT:
% beats: matrix with beat indexes, one set per column, filled with NaN
% QRSLIMIT: admited diference between marks of the same beat in samples (optional)
%           QRSLIMIT= 250 samples by default
%
%OUTPUT:
% ind: indexes of the initial column vectors for each set, corresponding to the sincronized beats
%      each line of ind corresponds to a beat and each collumn to a set
%      if a beat is absent of a set
%
% Rute Almeida OCT08
% Last update: 27JUL2011
if nargin<1
    messages.errors=[messages.errors {'Fatal error in sincronizebeats: not enough inputs.'}];
    warning(char(messages.errors(end)))
    messages.errors_desc=[messages.errors_desc 'Mandatory inputs not defined.'];
    messages.status=0;
    return
elseif nargin<2
    QRSLIMIT=250;
    messages.setup.sincronizabeats.QRSLIMIT=QRSLIMIT;
elseif nargin<3
    messages.warnings=[];
    if ~isnumeric(QRSLIMIT)
        QRSLIMIT=250;
        messages.setup.sincronizabeats.QRSLIMIT=QRSLIMIT;
        messages.warnings=[messages.warnings {['Invalid QRSLIMIT:' num2str(QRSLIMIT) 'used.']}];
    end
end
if ~isfield(messages,'setup'), messages.setup=[]; end
if ~isfield(messages,'errors'), messages.errors=[]; end
if ~isfield(messages,'errors_desc'), messages.errors_desc=[]; end
if ~isfield(messages,'warnings'), messages.warnings=[]; end
if isfield(messages,'status')
    if messages.status~=1
        messages.warnings=[messages.warnings {['Initial status=' num2str(messages.status) '. status changed to 1']}];
    end
end
messages.status=1;

try
    nsets=size(beats,2);
    A=NaN*ones(length(beats)*nsets,3);
    indA=1;
    for s=1:nsets
        %aux=beats(~isnan(beats(:,s)),s);
        aux=beats(:,s);
        A(indA:(indA+length(aux)-1),1)=aux;
        A(indA:(indA+length(aux)-1),2)=s*ones(size(aux));
        A(indA:(indA+length(aux)-1),3)=1:length(aux);
        indA=indA+length(aux);
    end
    %A(isnan (A(:,1)),:)=[]; %25MAY2011
    [aux,iA]=sort(A(:,1)); %#ok<ASGLU>
    A=A(iA,:);
    beat_rep=[];
    repeated_pos=[];
    ibeat=1;
    nbeats=0;
    ind=NaN*ones(fix(1.5*length(beats)),nsets);
    mat_aux=((1:nsets)'*ones(1,nsets+5))';
%     tic
    while (ibeat<=length(A))
        aux=A(ibeat,1)+QRSLIMIT;
        %[miibeat,iibeat]=max(find(A((ibeat):min(ibeat+nsets-1, length(A)),1)<aux));
        % check how many beats are within the tol of QRSLIMIT
        [miibeat]=find(A((ibeat):min(ibeat+2*nsets, length(A)),1)<aux,1,'last');
        beat=A(ibeat:(ibeat+miibeat-1),:);
        %if a set of marks is duplicated there are overlaping sets
        if miibeat>nsets %
            repeated_set=find(sum(beat(:,2)*ones(1,nsets)-mat_aux(1:miibeat,:)==0)>1);
            n_rep=size(repeated_set,2);
            for k=1:n_rep
                repeated_pos = [repeated_pos find(beat(:,2)== repeated_set(k))]; %#ok<AGROW>
                [a,beat_c(k)]= min(abs((beat(repeated_pos(:,k),1)-nanmedian(beat(setdiff(1:size(beat,1),repeated_pos(:,k)),1))))); %#ok<ASGLU,AGROW>
                repeated_pos(beat_c(k),k)=NaN;
            end
            %an overlaping sets at the begging constitutes is own annotation
            if any(repeated_pos==1)
                beat_rep=beat(repeated_pos(~isnan(repeated_pos)),:);
            elseif any(repeated_pos<size(beats,2))   
                kindx = find(repeated_pos<size(beats,2));
                for k=1:length(find(repeated_pos<size(beats,2)))                    
                    messages.warnings=[messages.warnings {['Mixed overlap: enable to sincronize correctely  beat ' num2str(nbeats+1) '. The beat ' num2str(beat(repeated_pos(kindx(k)),3)) ' of the set ' num2str(beat(repeated_pos(kindx(k)),2)) 'was excluded.']}];
                end
            end
            beat(repeated_pos(~isnan(repeated_pos)),:)=[];
            repeated_pos=[];
        end
        if ~isempty(beat_rep)
            if all(beat_rep(:,2)==unique(beat_rep(:,2))); %for overlaping sets at the begging, check for new overlaps
                nbeats=nbeats+1;
                ind(nbeats,beat_rep(:,2))= beat_rep(:,3);
            else
                for k=1:n_rep
                    nbeats=nbeats+1;
                    ind(nbeats,beat_rep(k,2))= beat_rep(k,3);
                    messages.warnings=[messages.warnings {['Doble overlap: enable to sincronize correctely beat' num2str(nbeats)]}];
                end
            end
        end
        nbeats=nbeats+1;
        ind(nbeats,beat(:,2))= beat(:,3);
         
        ibeatf=ibeat+size(beat,1)+size(beat_rep,1)-1;
        ibeat=max(ibeatf+1,ibeat+1);
        beat_rep=[];
        n_rep=0;
    end
%     toc
    ind(ind==0)=NaN;
    aux=find(diff(ind(:,1))==0);
    if( ~isempty( aux ) )
        ind(aux(isnan(ind(aux+1,2)))+1,1)=NaN;%FEB15_2012
    end
    ind(sum(isnan(ind),2)==size(beats,2),:)=[];
    
catch me
    messages.errors=[messages.errors {'Fatal error in sincronizebeats.'}];
%     warning(char(messages.errors(end)))
    messages.errors_desc=[messages.errors_desc {me.message}];
    messages.status=0;
    return
end
