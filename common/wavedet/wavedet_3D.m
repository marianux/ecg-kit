function [position,positionaux,messages] = wavedet_3D(sig_all, ext_anot, heasig, messages )
%[position,positionaux,messages] =wavedet_3D(sigdir,headir,matdir,ecgnr,ft,anot,lead,t,qrs_flag,anot_fmt,aname,dirann)
%Significant ECG points detector based on wavelet approach
%
%Input Parameters:
%	sigdir: directory of the original signal (.dat file)
%	headir: directory of the header file
%	matdir: output directory
%	ecgnr: string with the name of the record to analyze
% 	ft: format file:
%           - 0 for MIT header using any recorded leads from same datafile
%           - 1 for LUND header (to be tested)
%           - 2 for data in matlab file (.mat file with a varible fa with sampling frequency and variable sinal with ecg, one lead per column)
%           - 4 for RawData.bin files (one hour) created with Mortara Software SuperECG
%
%          obsolete formats
%           - 3 for MIT header using the Dower xyz (ft=0, leadsynth_flag=1)
%           - 10 for MIT format in PTBDB .dat extension file (ft=0)
%           - 30 for MIT type format in PTBDB .dat extension file use the projected xyz (ft=0,  leadsynth_flag=1)
%           - 20 for MIT type format in PTBDB .xyz extension file (ft=0, leads=leads+12)
%           - 40 for MIT header using 3 PC based all available leads and signal (ft=0, leadsynth_flag=4,flagaux=0)
%           - 41 for MIT header using 3 PC based in 8 leads and sub interval on the beat (ft=40,leadsynth_flag=3,flagaux=1)
%           - 42 for MIT type format in PTBDB .dat extension file use the PC based all signal and 8 leads (ft=0,leadsynth_flag=3,flagaux=0)
%           - 44 for MIT header using for MIT type using lead aVF and 2 PC out of precordial leads and all signal (ft=0,leadsynth_flag=5,flagaux=0)
%           - 123 for RawData.bin files (one hour) created with Mortara Software SuperECG (ft=4)
%           - 124 for RawData.bin files (one hour) created with Mortara Software SuperECG use the projected xyz (ft=4, leadsynth_flag=1)
%           - 50 for MIT type format like in ft=0 file using the PC based all available leads and signal (ft=0, leadsynth_flag=4)
%           - 51 for MIT type format like in ft=0 file using the  PC based in 8 leads and sub interval on the beat (ft=0,leadsynth_flag=3,flagaux=1)
%           - 52 for MIT type format like in ft=0 file use the PC based all signal and 8 leads (ft=40,leadsynth_flag=3,flagaux=0)
%
%	anot: name of the annotation output file
%	lead: vector with the leads number to analyze (1,2 or 3 leads)
%         if lead is a cell array the first element stands for leads number to analyze from the ones read and
%         the second for the leads to be read from the original signal file
%	t=[tbegin tend]: time vector with initial and sample to analyze
%   flags: [qrs_flag leadsynth_flag flagaux]
%           qrs_flag: QRS detection only (2), External (1), internal (0) QRS detector (0 by default)
%           leadsynth_flag: None (0), (default)
%                           Dower matrix derived VCG (1)
%                           Levkov's matrix derived VCG (2)
%                           Principal components based in 8 uncorrelated leads out of 12-lead standard (3)
%                           Principal components based in all available leads (4)
%                           lead aVF and Principal components based on precordial leads (5)
%           flagaux: use all signal (0) (default) or sub interval on the beat (1) for PC calculation, available only for leadsynth_flag>2
%
%   anot_fmt: cell with format in case of external annotation file
%             first field for QRS detector, format of annotation file: MIT (0) or mat file (1)
%             second field for with excluded segments: xls (0), MIT annotations (1)
%   aname: cell with annotation file name with QRS marks in the first field (by default ecgnr)
%               with annotation file name with excluded segments in the second field (by default aname{1})
%   dirann: cell with directory with external annotation files (QRS marks and/or excluded segments, by default sigdir)
%          if diferent annotation files are used, first field refers to QRS marks and second field to excluded segments
%
%Output Parameters:
%	position: struct vector with the detected points locations in samples
%
%   positionaux: struct vector with the detected points using each specified lead
%                   - positionaux.position1: using lead(1)
%                   - positionaux.position2: using lead(2)
%                   - positionaux.position3: using lead(3)
%
%  messages: struct vector with errors, warnings, setup and state
%                   - messages.status=0,1;
%                   - messages.errors
%                   - messages.errors_desc
%                   - messages.warnings
%                   - messages.setup.wavedet.nsamp=nsamp; number of samples to process in each excerpt, consider to calculate some parameters (default 2^16/250*sf)
%                   - messages.setup.wavedet.sig_quality: sig_quality test 1 or 0 (activated or inactivated)
%
% struct description:
%             Pon: P wave onset
%              P: P wave peak
%           Poff: P wave end
%          QRSon: QRS complex onset
%              Q: Q wave peak (in multilead, according to QRSonset best fitted lead)
%      R_inQRSon: R wave peak in QRSonset best fitted lead (multilead only)
%              R: R wave peak (median mark in multilead approach)
%            qrs: QRS complex fiducial mark
%     R_inQRSoff: R wave peak in QRSend best fitted lead (multilead only)
%         Rprima: R' wave peak (in multilead, according to QRSend best fitted lead)
%              S: S wave peak (in multilead, according to QRSend best fitted lead)
%         QRSoff: QRS complex end
%            Ton: T wave onset
%              T: first T wave peak (median mark in multilead approach)
%         Tprima: second T wave peak (in multilead, according to Tend best fitted lead)
%           Toff: T wave end
%          Ttipo: T wave morphology (1:5 corresponding respectively to normal, inverted, upwards only, downwards only, biphasc pos-neg, biphasc neg-pos)
%         Tscale: scale used for T wave detection
%        Ttipoon: T wave morphology, according to Tonset best fitted lead (multilead only)
%       Ttipooff: T wave morphology, according to Tend best fitted lead (multilead only)
%
% warnings description:
%       'warning: there are not 8 uncorrelated leads, all available leads considered instead'
%       leadsynth_flag=3 was used over a file that does not have all precordial leads and at least 2 of the limb leads
%       leadsynth_flag=4 was automatically set
%
%       'warning: lead aVF not available, lead xpto considered instead'
%       leadsynth_flag=5 was used over a file that does not have a lead aVF nor 2 limb leads to calculate it
%       aVF was replaced by lead xpto
%
%       'warning: mV are assumed as units'
%       not units are provided, mV are assumed
%
% The analysis loop has an overlapping structure.  The number of samples
% processed is large and there are an overlap at the beginning and another
% at the end.  The first one is necessary because: a) the first l4-1 samples
% are always incorrectly filtered b) To align the signals, as the filters
% have different lengths, we must discard some samples. c) The algorithms
% must have some possibility of turning back.
% The second overlap is mainly because of the alignment, as the filters have
% different delays.
%
% based on wavedet.m by Juan Pablo Mart�nez Cort�s
% previous versions include wavedetplus.m and wavedet3D.m
% Last update: Rute Almeida  07FEB2012
%
% Designed for MATLAB Version R12; tested with MATLAB Version R13

ft = 2;
str = '1';
positionaux=[];
position=[];
aname_sig_quality_cell{3} = [];
t_initial=0;  % Juan 28/03/2012
lead = 1;

% if (nargin<14) % JB 27/07/2011
%     %     bufferedSignal = [];  % Juan 28/03/2012  No se usan no tienen porque
%     %     estar definidas
%     %     bufferedHeaSig = [];  % Juan 28/03/2012
% else % 03AGO 2011
%     ft=NaN;
%     heasig = bufferedHeaSig;
% %     if ~isfield(heasig,'freq')
% %         messages.errors=[messages.errors {'Fatal error in wavedet_3D: no sampling frequency variable found.'}];
% %         warning(char(messages.errors(end)))
% %         messages.errors_desc=[messages.errors_desc {'No sampling frequency variable found on the header information file.'}];
% %         messages.status=0;
% %         return
% %     end
% %     if ~isempty(t) && length(t)==2 %& t(2)>t(1)
% %         if length(bufferedSignal) == t(2)-t(1)+1;
% %             sig_all =bufferedSignal(:,~isnan(bufferedSignal(1,:)));
% %             %             t_initial=t(1);  % Juan 28/03/2012
% %             t=[1 length(sig_all)];
% %         elseif length(bufferedSignal)>t(2)
% %             sig_all = bufferedSignal(:,t(1):t(2));
% %         else
% %             sig_all = bufferedSignal;
% %             t=[1 length(sig_all)];
% %         end
% %     else
% %         sig_all = bufferedSignal;
% %         t=[1 length(sig_all)];
% %     end
% %     
%     
%     
%     if size(sig_all,2)>size(sig_all,1)
%         sig_all=sig_all';
%     end
%     
%     if ~isfield(heasig,'nsamp') || isempty(heasig.nsamp)
%         heasig.nsamp=size(sig_all,1);
%     end
%     if ~isfield(heasig,'nsig')
%         heasig.nsig=size(sig_all,2);
%     end
%     if ~isfield(heasig,'gain')
%         heasig.gain=200.*ones(1,heasig.nsig);  % When heasig.gain = 0 => default 200
%     end
%     if ~isfield(heasig,'spf')
%         %         heasig.spf=ones(1,max(heasig.nsig,leadreading));  % Juan
%         %         28/03/2012  leadreading no esta definido
%         heasig.spf=ones(1,max(heasig.nsig, size(sig_all,1)));
%     end
%     if ~isfield(heasig,'adczero')
%         heasig.adczero=zeros(1,heasig.nsig);
%     end
% end

t=[1 size(sig_all,1)];

%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    messages.warnings=[];
end

if ~isfield(messages,'setup'), messages.setup.wavedet=[]; end
if ~isfield(messages.setup.wavedet,'QRS_detection_only'), messages.setup.wavedet.QRS_detection_only=false; end
if ~isfield(messages.setup.wavedet,'QRS_detection_thr'), messages.setup.wavedet.QRS_detection_thr=ones(5,1); end
if ~isfield(messages,'errors'), messages.errors=[]; end
if ~isfield(messages,'errors_desc'), messages.errors_desc=[]; end
if ~isfield(messages,'warnings'), messages.warnings=[]; end
if isfield(messages,'status')
    if messages.status~=1
        messages.warnings=[messages.warnings {['Initial status=' num2str(messages.status) '. status changed to 1']}];
    end
end
messages.status=1;
messages.setup.wavedet.sig_quality = 2;

% if nargin<7 %21MAR2011
%     messages.errors=[messages.errors {'Fatal error in wavedet_3D: not enough inputs.'}];
%     warning(char(messages.errors(end)))
%     messages.errors_desc=[messages.errors_desc 'Mandatory inputs not defined.'];
%     messages.status=0;
%     return
% elseif nargin<9
%     flags=[];
% end
% if iscell(lead)
%     if size(lead,2)==2
%         leadreading=(lead{2});
%     else
%         leadreading=(lead{1});
%     end
%     lead=lead{1};
% else
%     leadreading=(lead);
% end

% if nargin<=10 %&& anot_fmt==1 %JUL2011
%     aname=ecgnr;
%     dirann=sigdir;
% end%JUL2011
% if iscell(aname)
%     if length(aname)==2
%         aname_sig_quality=aname{2};
%     end
%     aname=aname{1};
% else
%     aname_sig_quality=aname;
% end
% if iscell(dirann)
%     if length(dirann)==2
%         dirann_sig_quality=dirann{2};
%     end
%     dirann=dirann{1};
% else
%     dirann_sig_quality=dirann;
% end

% if isempty(flags)
%     if isfield(messages.setup.wavedet,'qrs_flag')
%         flags(1)=messages.setup.wavedet.qrs_flag;
%     else
%         flags(1)=0;
%     end
% elseif isfield(messages.setup.wavedet,'qrs_flag') &&  flags(1)~=messages.setup.wavedet.qrs_flag;
%     messages.warnings=[messages.warnings {'Obsolete use of flag with inconsistency: flag(1) value differs from messages.setup.wavedet.qrs_flag value. flag(1) ignored.'}];
% else
%     messages.setup.wavedet.qrs_flag=flags(1);
%     messages.warnings=[messages.warnings {['Obsolete use of flag. Use messages.setup.wavedet.qrs_flag=' num2str(flags(1)) 'instead.' ]}];
% end
% if length(flags)<2
%     if isfield(messages.setup.wavedet,'leadsynth_flag')
%         flags(2)=messages.setup.wavedet.leadsynth_flag;
%     else
%         flags(2)=0;
%     end
% elseif isfield(messages.setup.wavedet,'leadsynth_flag') &&  flags(2)~=messages.setup.wavedet.leadsynth_flag;
%     messages.warnings=[messages.warnings {'Obsolete use of flag with inconsistency: flag(2) value differs from messages.setup.wavedet.leadsynth_flag value. flag(2) ignored.'}];
% else
%     messages.setup.wavedet.leadsynth_flag=flags(2);
%     messages.warnings=[messages.warnings {['Obsolete use of flag. Use messages.setup.wavedet.leadsynth_flag=' num2str(flags(2)) 'instead.' ]}];
% end
% if length(flags)<3
%     if isfield(messages.setup.wavedet,'flagaux')
%         flags(3)=messages.setup.wavedet.flagaux;
%     else
%         flags(3)=0;
%     end
% elseif isfield(messages.setup.wavedet,'flagaux') &&  flags(3)~=messages.setup.wavedet.flagaux;
%     messages.warnings=[messages.warnings {'Obsolete use of flag with inconsistency: flag(3) value differs from messages.setup.wavedet.flagaux. flag(3) ignored.'}];
% else
%     messages.setup.wavedet.flagaux=flags(3);
%     messages.warnings=[messages.warnings {['Obsolete use of flag. Use messages.setup.wavedet.flagaux=' num2str(flags(3)) 'instead.' ]}];
% end
% if length(flags)==4
%     if isfield(messages.setup.wavedet,'nsamp') && messages.setup.wavedet.nsamp~=flags(4)
%         messages.warnings=[messages.warnings {'Obsolete use of flag with inconsistency: flag(4) value differs from messages.setup.wavedet.nsamp value. flag(4) ignored.'}];
%     else
%         messages.setup.wavedet.nsamp=flags(4);
%         messages.warnings=[messages.warnings {['Obsolete use of flag. Use messages.setup.wavedet.nsamp=' num2str(flags(4)) 'instead.' ]}];
%     end
% end


if( nargin < 2 || isempty(ext_anot) )
    qrs_flag=0;
else
    qrs_flag=1;
end

leadsynth_flag=0;
flagaux=0;

flags = [ qrs_flag leadsynth_flag flagaux];

messages.setup.wavedet.qrs_flag=qrs_flag;
messages.setup.wavedet.leadsynth_flag=leadsynth_flag;
messages.setup.wavedet.flagaux=flagaux;

% if ~isfield(messages.setup.wavedet,'sig_quality') || messages.setup.wavedet.sig_quality==0
%     messages.setup.wavedet.sig_quality=0;
%     messages.warnings=[messages.warnings {'No signal quality test used: messages.setup.wavedet.sig_quality=0'}];
% else
%     if isfield(messages.setup.wavedet,'dirann_sig_quality')
%         dirann_sig_quality=messages.setup.wavedet.dirann_sig_quality;
%     elseif ~exist('dirann_sig_quality','var') || isempty(dirann_sig_quality)
%         dirann_sig_quality=sigdir;
%     end
%     if isfield(messages.setup.wavedet,'aname_sig_quality')
%         aname_sig_quality=messages.setup.wavedet.aname_sig_quality;
%     elseif ~exist('aname_sig_quality','var') || isempty(aname_sig_quality)
%         aname_sig_quality=ecgnr;
%     end
%     messages.setup.wavedet.sig_quality=1;
%     messages.warnings=[messages.warnings {['Signal quality information used: messages.setup.wavedet.sig_quality='  num2str(messages.setup.wavedet.sig_quality)]}];
% end

if length(flags)==3
    if  leadsynth_flag<3
        flagaux=0;
        messages.warnings=[messages.warnings {'flagaux>0 only available for leadsynth_flag>2. flagaux=0 used'}];
    end
end

%obsolete formats
% if ft==3;
%     ft=0; leadsynth_flag=1;
% end %obsolete formats
% if ft==10
%     ft=0;
% end %obsolete formats
% if ft==20
%     ft=0;
%     lead=lead+12;
%     leadreading=leadreading+12;
% end %obsolete format
% if ft==30
%     ft=0;
%     leadsynth_flag=1;
%     leadreading=1:12;
% end%obsolete formats
% if ft==41
%     ft=0;
%     leadsynth_flag=3;
%     flagaux=1;
% end %obsolete format
% if ft==42||ft==52
%     ft=0;
%     leadsynth_flag=3;
%     flagaux=0;
% end%obsolete formats
% if ft==50
%     ft=0;
%     leadsynth_flag=4;
%     flagaux=0;
% end%obsolete formats
% if ft==51
%     ft=0;
%     leadsynth_flag=4;
%     flagaux=1;
% end%obsolete formats
% if ft==123
%     ft=4;
% end% change format notation
% if ft==124
%     ft=4;
%     leadsynth_flag=1;
% end% change format notation
% if ft==44
%     ft=0;
%     leadsynth_flag=5;
%     flagaux=0;
% end%obsolete formats
% if ft==40
%     ft=0;
%     leadsynth_flag=4;
%     flagaux=0;
% end%obsolete formats
% 
% %number of leads from 1 to 3
% if length(lead)>3
%     lead(4:end)=[];
%     messages.warnings=[messages.warnings {['Too many leads to process, just first 3 leads are considered: lead=[' num2str(lead) ']']}];
% elseif isempty(lead)
%     lead=1;
%     leadreading=lead;
%     messages.warnings=[messages.warnings {'No leads indicated to process: first lead in file is considered'}];
% end
% 
% if leadsynth_flag==1
%     leadstring='synt';
% end
% if leadsynth_flag==2
%     leadstring='levk';
% end
% if leadsynth_flag==1 || leadsynth_flag==2
%     %%%%%%%%%%%%%%%%%%%%%%%
%     %verify leads to process
%     invalid_lead= find(lead>3); %#ok<EFIND>
%     if ~isempty(invalid_lead)
%         [C,IA,IB]=intersect(leadreading,lead); %#ok<NASGU> % Juan 28/03/2012
%         if length(C)==length(lead)
%             lead=IA;
%             invalid_lead= find(lead>3); % Juan 28/03/2012
%             if isempty(invalid_lead)
%                 messages.warnings=[messages.warnings {['the lead(s) chosen to be processed are the lead(s) ' num2str(leadreading(lead))]}];
%             else
%                 alead=lead;
%                 lead(invalid_lead)=[];
%                 if isempty(lead)
%                     lead=1:min(length(lead),3);
%                 end
%                 messages.warnings=[messages.warnings {[ 'Invalid lead for VCG analysis (' num2str(alead) '): transformed lead(s) ' num2str(1:min(length(lead),3)) ' considered instead.' ]}];
%             end
%         elseif ~isempty(C)
%             messages.warnings=[messages.warnings {[ 'Invalid lead for VCG analysis (' num2str(lead) '): transformed lead(s) ' num2str(leadreading(IA)) ' considered instead.' ]}];
%             lead=IA;
%         else
%             messages.warnings=[messages.warnings {[ 'Invalid lead for VCG analysis (' num2str(lead) '): transformed lead(s) ' num2str(1:min(length(lead),3)) ' considered instead.' ]}];
%             lead=1:min(length(lead),3);
%         end
%     end
%     %%%%%%%%%%%%%%%%
% elseif leadsynth_flag==3
%     invalid_lead=find(lead>8);
%     lead(invalid_lead)=[];
%     if ~isempty(invalid_lead)
%         messages.warnings=[messages.warnings {'Invalid lead after the chosen PC transformation: 8 PC are produced.'}];
%     end
% elseif leadsynth_flag==4
%     invalid_lead= find(lead>max(leadreading)); % Juan 28/03/2012
%     lead(invalid_lead)=[];
%     if ~isempty(invalid_lead)
%         messages.warnings=[messages.warnings {['Invalid lead after the chosen PC transformation: ' num2str(leadreading) ' PC are produced. Lead(s) ' num2str(lead)  ' considered instead']}];
%     end
%     if isempty(lead)
%         lead=1;
%         messages.warnings=[messages.warnings {['Invalid lead after the chosen PC transformation: lead=[' num2str(lead) '] considered instead.']}];
%     end
%     
% elseif leadsynth_flag==5
%     if isempty(find(lead==1,1))
%         messages.warnings=[messages.warnings { 'Processing lead 1, standing for aVF, is mandatory in the chosen PC transformation.'}];
%         lead=[1 lead];
%     end
%     invalid_lead=find(lead>7);
%     lead(invalid_lead)=[];
%     if ~isempty(invalid_lead)
%         messages.warnings=[messages.warnings {'Invalid lead after the chosen PC transformation: aVF plus 6 PC are available.'}];
%     end
%     if isempty(lead)
%         lead=1:3;
%         messages.warnings=[messages.warnings {['Invalid lead after the chosen PC transformation: lead=[' num2str(lead) '] considered instead.']}];
%     elseif length(lead)>3
%         lead(4:end)=[];
%         messages.warnings=[messages.warnings {['Too many leads to process, just first 3 leads are considered: lead=[' num2str(lead) ']']}];
%     elseif length(lead)==1
%         messages.warnings=[messages.warnings { 'Single lead over lead aVF selected.' }];
%     end
% end

% global OPT
allfig=0;figuresoff=1; %#ok<NASGU>

% OPT=optimset(@fminbnd);
% OPT.MaxFunEvals=1000;
% OPT.Display='off';
% warning off %#ok<WNOFF>
ultimo_anot=1;
ultimo_anottyp=0;
% timeoffset=0;   % Juan 28/03/2012  Solo se define para el formato lund si
% es necesario dentro
count=0;

%errores=struct('errores',[]);

sep = filesep;

% if isunix,          %%%%%%% Rute 13/03/02
%     sep = '/';
% else sep = '\';     %%%%%%% Rute 13/03/02
% end                 %%%%%%% Rute 13/03/02

PC_aux_ind=[];
% aux=find(ecgnr==sep);
% headir=[headir ecgnr(1:aux)];
% sigdir=[sigdir ecgnr(1:aux)];
% matdir=[matdir ecgnr(1:aux)];
% if ~isempty(aux)
%     ecgnr=ecgnr(aux+1:end);
% end

%%%%%%%%%%%%%%%%%%%%%%% file format %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if ft==0   % MIT format files
%     name=[headir sep ecgnr '.hea'];
%     heasig = readheader(name);
%     [auxleads] = ecgleads(heasig.desc); %V1 V2 V3 V4 V5 V6 I II III aVR aVL aVF X Y Z VX VY VZ
%     
%     %verify leads to read
%     if sum(leadreading>heasig.nsig)>0
%         messages.warnings=[messages.warnings {[ 'The selected file contains only ' num2str(heasig.nsig) ' signals: unable to read leads [' num2str(leadreading) '], leads out of the range will be ignored' ]}];
%         leadreading(leadreading>heasig.nsig)=[];
%     end
%     
%     if leadsynth_flag==0 % should be leadreading=lead
%         if length(leadreading)~=length(lead) | (length(leadreading)==length(lead) & sort(leadreading)~=sort(lead)) %#ok<AND2,OR2>
%             [C,IA,IB]=intersect(leadreading,lead); %#ok<NASGU>
%             if length(C)==length(lead) % extra leads on leadreading will be deleted
%                 leadreading=IA;
%                 messages.warnings=[messages.warnings {['the lead(s) chosen to be processed are the lead(s) ' num2str(leadreading)]}];
%             else % extra leads on lead or different leads
%                 if ~isempty(lead>heasig.nsig)
%                     messages.warnings=[messages.warnings {[ 'The selected file contains only ' num2str(heasig.nsig) ' signals: unable to read leads [' num2str(lead) '], leads out of the range will be ignored' ]}];
%                     lead=lead(lead<=heasig.nsig);
%                 end
%                 if length(leadreading)~=length(lead) | (length(leadreading)==length(lead) & sort(leadreading)~=sort(lead)) %#ok<AND2,OR2>
%                     messages.warnings=[messages.warnings {[ 'unable to interpret the lead(s) number to process (' num2str(leadreading) '): lead(s) ' num2str(lead(lead<=heasig.nsig)) ' considered instead.' ]}];
%                     leadreading=lead;
%                 end
%             end
%         else
%             messages.warnings=[messages.warnings {['the lead(s) chosen to be processed are the lead(s) ' num2str(leadreading)]}];
%         end
%     elseif leadsynth_flag==1 || leadsynth_flag==2 || leadsynth_flag==3% Leads V1 to V6 and 2 limb leads required
%         if size(leadreading)<8 %wrong number of leads
%             if ~isnan(auxleads(1:8))
%                 leadreading=(auxleads(1:8));
%             elseif ~isnan(auxleads([1:7 9]))
%                 leadreading=(auxleads([1:7 9]));
%             elseif ~isnan(auxleads([1:6 8 9]))
%                 leadreading=(auxleads([1:6 8 9]));
%             else
%                 messages.errors=[messages.errors {'Fatal error in wavedet_3D: not enough leads in file to apply the chosen transformation.'}];
%                 warning(char(messages.errors(end)))
%                 messages.errors_desc=[messages.errors_desc 'Leads V1 to V6 and 2 limb leads are required.'];
%                 messages.status=0;
%                 return
%             end
%             messages.warnings=[messages.warnings {['Not enough leads to apply the chosen transformation in given leadreading: leadreading=['  num2str(leadreading) '] considered instead' ]}];
%         else
%             if sum(isnan(auxleads(1:6)))>0 || sum(isnan(auxleads(7:9)))>1 % heasig.desc does not include Leads V1 to V6 and 2 limb leads
%                 if heasig.nsig>=8
%                     messages.warnings=[messages.warnings {'Leads V1 to V6 and 2 limb leads required to apply the chosen transformation but cannot be recognized. It is assumed that the first 8 leads of leadreading are V1 V2 V3 V4 V5 V6 I II'}];
%                     leadreading=leadreading(1:8);
%                     auxleads=NaN*ones(size(auxleads));
%                     auxleads(1:8)=1:8;
%                 else
%                     messages.errors=[messages.errors {'Fatal error in wavedet_3D: not enough leads in file to apply the chosen transformation.'}];
%                     warning(char(messages.errors(end)))
%                     messages.errors_desc=[messages.errors_desc 'Leads V1 to V6 and 2 limb leads are required.'];
%                     messages.status=0;
%                     return
%                 end
%             else % heasig.desc includes Leads V1 to V6 and 2 limb leads
%                 leadreading=leadreading(auxleads(auxleads<=length(leadreading)));% sort leadreding as V1 V2 V3 V4 V5 V6 I II III aVR aVL aVF X Y Z VX VY VZ
%             end
%         end
%     elseif leadsynth_flag==4 && size(leadreading)==1
%         leadsynth_flag=0;
%         messages.warnings=[messages.warnings {'To aply PC transformation more than one lead is requires in leadreading: no trasnformation applied (leadsynth_flag=0)'}];
%     elseif leadsynth_flag==5 && length(lead)~=1 %  aVL and leads V1 to V6 required
%         if size(leadreading)<7 %wrong number of leads
%             if ~isnan(auxleads([12 1:6])) %lead aVF and V1 to V6 available
%                 leadreading=auxleads([1:6 12]);
%             elseif ~isnan(auxleads(1:8))
%                 leadreading=auxleads(1:8);
%             elseif ~isnan(auxleads([1:7 9]))
%                 leadreading=auxleads([1:7 9]);
%             elseif ~isnan(auxleads([1:6 8 9]))
%                 leadreading=auxleads([1:6 8 9]);
%             else
%                 messages.errors=[messages.errors {'Fatal error in wavedet_3D: not enough leads in file to apply the chosen transformation.'}];
%                 warning(char(messages.errors(end)))
%                 messages.errors_desc=[messages.errors_desc 'Leads aVF and V1 to V6 are required. Sugestion: use leadsynth_flag=4 instead'];
%                 messages.status=0;
%                 return
%             end
%         else
%             if sum(isnan(auxleads(1:6)))>0 || (isnan(auxleads(12)) && sum(isnan(auxleads(7:9)))>1) % heasig.desc does not include Leads V1 to V6 and aVF cannot be calculated
%                 if size(leadreading)==7 &&  heasig.nsig>=7
%                     messages.warnings=[messages.warnings {'Leads V1 to V6 and aVF limb leads required to apply the chosen transformation but cannot be recognized. It is assumed that the 7 leads of leadreading are V1 V2 V3 V4 V5 V6 aVF' }];
%                     leadreading=leadreading(1:7);
%                     auxleads=NaN*ones(size(auxleads));
%                     auxleads([1:6 12])=1:7;
%                 elseif heasig.nsig>=8
%                     messages.warnings=[messages.warnings {'Leads V1 to V6 and aVF limb leads required to apply the chosen transformation but cannot be recognized. It is assumed that the first 8 leads of leadreading are V1 V2 V3 V4 V5 V6 aVF' }];
%                     leadreading=leadreading(1:8);
%                     auxleads=NaN*ones(size(auxleads));
%                     auxleads(1:8)=1:8;
%                 else
%                     messages.errors=[messages.errors {'Fatal error in wavedet_3D: not enough leads in file to apply the chosen transformation.'}];
%                     warning(char(messages.errors(end)))
%                     messages.errors_desc=[messages.errors_desc 'Leads aVF and V1 to V6 are required. Sugestion: use leadsynth_flag=4 instead'];
%                     messages.status=0;
%                     return
%                 end
%             else % heasig.desc includes Leads V1 to V6 and 2 limb leads
%                 leadreading=leadreading(auxleads(auxleads<=length(leadreading)));% sort leadreding as %V1 V2 V3 V4 V5 V6 I II III aVR aVL aVF X Y Z VX VY VZ
%             end
%         end
%     elseif leadsynth_flag==5 && length(lead)==1
%         if ~isnan(auxleads(12))
%             leadreading=auxleads(12);
%         elseif sum(~isnan(auxleads(7:9)))>1
%             leadreading=auxleads(~isnan(auxleads(7:9)));
%         elseif length(leadreading)==1
%             messages.warnings=[messages.warnings {'Lead aVF lead is required but cannot be recognized. It is assumed that the first lead of leadreading is aVF'}];
%             leadreading=1;
%             auxleads=NaN*ones(size(auxleads));
%             auxleads(12)=1;
%         else
%             messages.warnings=[messages.warnings {'Lead aVF lead is required but cannot be recognized. It is assumed that the first 2 lead of leadreading are I and II'}];
%             auxleads=NaN*ones(size(auxleads));
%             auxleads(7:8)=1:2;
%         end
%     end
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% this is necessary when there are leads in different files.   e.g. ptbdb database  JP 2006
%     leadsfile=NaN*ones(heasig.nsig,1);
%     leadinfile=NaN*ones(heasig.nsig,1);
%     nsiginfile=NaN*ones(heasig.nsig,1);
%     for jj=1:heasig.nsig
%         for kk=1:heasig.nsig
%             leadsfile(kk)=strcmp(heasig.fname(jj,:),heasig.fname(kk,:));
%         end
%         leadsfile = cumsum(leadsfile);
%         leadinfile(jj) = leadsfile(jj);  % lead indexes in the file
%         nsiginfile(jj) = leadsfile(end); % number of the leads in th file nsiginfile<=heasig.nsig
%     end
%     
%     if(min(lead)>nsiginfile(1))
%         auxleads(auxleads<=nsiginfile(1))=NaN;
%     else
%         auxleads(auxleads>nsiginfile(1))=NaN;
%     end
%     n_auxleads=find(auxleads>0);
%     auxleads(n_auxleads)=leadinfile(auxleads(n_auxleads));
%     %     nsiginfile=nsiginfile(lead(1)); % Juan 28/03/2012
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%
%     if (heasig.fmt(lead)==16) | (heasig.fmt(lead)==212) | (heasig.fmt(lead)==61), %#ok<OR2>
%         %%%%%%% multilead%%% Rute 02.Dec.04
%         if size(lead,2)==1
%             formato = num2str(heasig.fmt(lead));
%         elseif   diff(heasig.fmt(lead))==0
%             formato = num2str(heasig.fmt(lead(1)));
%         end
%         %%%%%%%%%%%%%
%     else
%         messages.errors=[messages.errors {'This format is not supported by the program.'}];
%         warning(char(messages.errors(end)))
%         messages.errors_desc=[messages.errors_desc ['Format ' num2str(heasig.fmt(lead)) ' data. For MIT data, formats 16, 212 and 61 are supported.'] ];
%         messages.status=0;
%         return
%     end
%     if strcmp(formato,'16')||strcmp(formato,'61'),
%         fid = fopen([sigdir heasig.fname(lead(1),:)],'rb');
%         fseek(fid,0,-1);  % Rewind the file
%         if strcmp(heasig.fname(1,1),'_'),  % Siemens card recordings with MIT-type header
%             timeoffset=512;
%             fseek (fid, timeoffset,-1); % offset
%         end
%     end
%     
% elseif ft==1  % LUND format files
%     try
%         hdsig=gethdsig(sigdir,ecgnr); % Reading header information
%         heasig = hdsig2heasig(hdsig);
%         leadinfile=lead;
%         
%         [auxleads] = ecgleads(heasig.desc);
%     catch me
%         messages.errors=[messages.errors {'Fatal error in wavedet_3D: unable to read header signal on Lund format.'}];
%         warning(char(messages.errors(end)))
%         messages.errors_desc=[messages.errors_desc {me.message}]; % Juan 28/03/2012
%         messages.status=0;
%         return
%     end
% elseif ft==2  % data in a Matlab file
%     try
%         [sig_all,messages1,heasig] = readsignal(sigdir,ecgnr,t,ft,[],leadreading,messages);  % Juan 28/03/2012
%         if messages1.status==0
%             error('Fatal error in wavedet_3D: unable to read signal on mat format.');  % Juan 28/03/2012
%         end
%     catch me
%         messages.errors=[messages.errors {messages1.errors}];  % Juan 28/03/2012
%         warning(char(messages.errors(end)))
%         messages.errors_desc=[messages.errors_desc {me.message}];  % Juan 28/03/2012
%         messages.status=0;
%         return
%     end
%     if size(sig_all,2)>size(sig_all,1)
%         sig_all=sig_all';
%     end
%     %should not be needed;   I will remove that part, all of that is done in
%     %readsignal  Juan 28/03/2012
%     if ~exist('heasig','var')
%         heasig.nsamp=size(sig_all,1);
%         try
%             if exist('fa','var')
%                 heasig.freq=fa; %#ok<NODEF>
%             else
%                 heasig.freq=messages.setup.wavedet.freq;
%             end
%         catch me
%             messages.errors=[messages.errors {'Fatal error in wavedet_3D: no sampling frequency variable found.'}];
%             warning(char(messages.errors(end)))
%             me.message = 'No sampling frequency variable found on the mat file.';  % Juan 28/03/2012
%             messages.errors_desc=[messages.errors_desc {me.message}];  % Juan 28/03/2012
%             messages.status=0;
%             return
%         end
%         heasig.nsig=size(sig_all,2);
%         heasig.gain=200.*ones(1,heasig.nsig);  % When heasig.gain = 0 => default 200
%         heasig.adczero=zeros(1,heasig.nsig);
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%
%     
%     leadinfile=lead;
%     %     nsiginfile=heasig.nsig;  %Juan 28/03/2012
%     heasig.spf_ecg=1;
%     
%     if isfield(heasig,'desc')
%         [auxleads] = ecgleads(heasig.desc);
%     else
%         auxleads = 1:heasig.nsig; % assuming that all are ECG
%     end
%     
% elseif ft==4  % Ana e Tiago 11-04-2006 (copy from wavedetplus)
%     fa=1000;
%     %     aux=[sigdir ecgnr '.bin'];  Juan 28/03/2012
%     if nargin >=8,
%         if isempty(t)
%             t=[1 fa*60*60];
%             %         else
%             %             t=[t(1) min(t(1)+fa*60*60,t(2))];
%         end
%     else
%         t=[1 fa*60*60];
%     end
%     %     first_sample=t(1); %ver com Sonia se '1 -1 pq  first_sample=0 por
%     %     defeito!!!!  Juan 28/03/2012
%     %     n_sample_per_lead=t(end);Juan 28/03/2012
%     %     lead_number=lead;Juan 28/03/2012
%     heasig.nsamp=t(2)-t(1);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     heasig.freq=fa;
%     heasig.gain=200;
%     heasig.nsig=12;
%     heasig.units=char(ones(heasig.nsig,1)*'microV');
%     heasig.desc=['I  ';'II ';'V1 ';'V2 ';'V3 ';'V4 ';'V5 ';'V6 ';'III';'aVR';'aVL';'aVF'];
%     [auxleads] = ecgleads(heasig.desc);
%     leadinfile=lead;
%     %     nsiginfile=heasig.nsig;
%     
% elseif isnan(ft) % data in buffer % 03AGO 2011
%     leadinfile=lead;
%     %     nsiginfile=heasig.nsig;           Juan 28/03/2012
%     heasig.spf_ecg=1;
%     if ~isfield(heasig,'nsamp') || heasig.nsamp<=0
%         heasig.nsamp=size(sig_all,1);
%     end
%     if isfield(heasig,'desc')
%         [auxleads] = ecgleads(heasig.desc);
%     else
%         auxleads = 1:heasig.nsig; % assuming that all are ECG
%     end
% else
%     try
%         [sig_all,messages1,heasig] = readsignal(sigdir,ecgnr,t,ft,[],leadreading,messages);
%         if messages1.status==0
%             error('Fatal error in wavedet_3D: no admissible format.');  % Juan 28/03/2012
%         end
%     catch me
%         messages.errors=[messages.errors {messages1.errors}];
%         me.message = 'Format ft should be one of {0,1,2,4}. See help wavedet3D for details.'; % Juan 28/03/2012
%         warning(char(messages.errors(end)))
%         messages.errors_desc=[messages.errors_desc {me.message}];  % Juan 28/03/2012
%         messages.status=0;
%         return
%     end
%     if size(sig_all,2)>size(sig_all,1)
%         sig_all=sig_all';
%     end
% end
% leadinfile=lead;  % Juan 28/03/2012  definido dentro de cada formato
% nsiginfile=heasig.nsig;  % Juan 28/03/2012 parece que no se usa m�s no es
% necesario definir
if ~isfield(heasig,'spf')
    heasig.spf=ones(1,heasig.nsig);
end
% try
%     heasig.spf_ecg=heasig.spf(leadinfile(1));
% catch me
%     me.message = 'It is taken the first lead in file.';
%     messages.warnings = [messages.warnings {me.message}];
%     heasig.spf_ecg=heasig.spf(1);
% end
heasig.spf_ecg = 1;
messages.setup.wavedet.freq=heasig.spf_ecg*heasig.freq; %13ENE2010
% if isfield(heasig,'desc')
%     [auxleads] = ecgleads(heasig.desc);
% else
%     auxleads = 1:heasig.nsig; % assuming that all are ECG
% end
% at this point heasig must exist
if isfield(messages.setup.wavedet,'nsamp')
    nsamp = messages.setup.wavedet.nsamp; % Number of samples per excerpt for the WT as input
else
    nsamp=0;
end
if nsamp==0
    nsamp = 2^16/250*messages.setup.wavedet.freq; % Number of samples per excerpt for the WT corresponding to 2^16 samples at sf=250
end
messages.setup.wavedet.nsamp=nsamp;
% if nargin >=8,
%     if isempty(t)
%         t=[1 heasig.nsamp]; % pode falhar em formato lund...
%     else
%         t(1)=max(t(1),1);
%         if ~isempty(heasig.nsamp)
%             t(2)=min(t(1)+heasig.nsamp-1,t(2));
%         end
%     end
% else
%     t=[1 heasig.nsamp];
% end
t=[1 heasig.nsamp];

% if flagaux==1
%     if isfield(messages.setup.wavedet,'leadsynth')
%         leadref= messages.setup.wavedet.leadsynth;
%         messages.warnings=[messages.warnings {'Lead' num2str(leadref) ' will be considered for PCA time interval.'}];
%     else
%         if ~isnan(auxleads(8)); %lead II
%             leadref=auxleads(8);
%         else
%             messages.warnings=[messages.warnings {'Lead II not available. First selected lead would be considered instead.'}];
%             warning(char(messages.warnings(end)));
%             leadref=lead(1);
%         end
%     end
%     if ~exist([matdir ecgnr '.s' num2str(leadref)],'file')
%         wavedet_3D(sigdir,headir,matdir,ecgnr,ft,['s' num2str(leadref)],leadref);
%     end
%     marks=readannot([matdir ecgnr '.s' num2str(leadref)]);
%     if exist('marks','var')
%         marks=posmat2(marks);
%         marks(marks==0)=NaN;
%     end
%     PC_aux_ind= [];
%     tol1=0;tol2=0; tol3=200;
%     
%     mark_on=4;mark_off=10;mark_on2=5; mark_off2=6;   % QRSonset-tend
%     %         mark_on=7;mark_off=10;mark_on2=6;mark_off2=10; % ton-tend
%     %         mark_on=8;mark_off=10;mark_on2=6;mark_off2=10;% tpeak-tend
%     %         mark_on=4;mark_off=6;mark_on2=4;mark_off2=6;% QRS based
%     
%     for tt=1:size(marks,1)
%         if ~isnan(marks(tt,mark_off)) & ~isnan(marks(tt,mark_on))%#ok<AND2> %
%             PC_aux_ind=[PC_aux_ind round((marks(tt,mark_on)-tol1)):round((marks(tt,mark_off)+tol2))]; %#ok<AGROW>
%         elseif ~isnan(marks(tt,mark_on))& ~isnan(marks(tt,mark_off2))%#ok<AND2> %
%             PC_aux_ind=[PC_aux_ind round((marks(tt,mark_on)-tol1)):(marks(tt,mark_off2)+tol3)]; %#ok<AGROW>
%         elseif ~isnan(marks(tt,mark_on2))& ~isnan(marks(tt,mark_off))%#ok<AND2> %
%             PC_aux_ind=[PC_aux_ind round((marks(tt,mark_on2)-tol3)):(marks(tt,mark_off)+tol2)]; %#ok<AGROW>
%         else
%             PC_aux_ind=[PC_aux_ind round((marks(tt,5)-tol3)):(marks(tt,5)+tol3)]; %#ok<AGROW>
%         end
%     end
% else
%     PC_aux_indaux= 1:(t(2)-t(1)+1);   %EXARLE OJO %PC_aux_indaux=[];%alternativa
% end
% try
%     heasig.gain(heasig.gain==0)=200;  % When heasig.gain = 0 => default 200
% catch me
%     me.message = 'Signal gain is set to 200'; % Juan 28/03/2012
%     messages.warnings = [messages.warnings {me.message}]; % Juan 28/03/2012
% end


% Reading of external QRS annotation file if qrs_flag
% if qrs_flag==1     % The QRS fiducial point is read from external annotator
%     switch (anot_fmt);
%         case 0        % MIT annotation file
%             if exist([dirann filesep aname],'file')%Jul2011 % Jbolea 01/12/11
%                 s=readannot([dirann aname],t);
%                 s=isqrs(s,heasig,t);
%                 ext_anot=s.time';
%             else error('QRS annotation file not found');
%             end
%         case 1        % mat file
%             if exist([dirann aname],'file')%Jul2011
% %                 if exist([dirann aname '.mat'],'file')%Jul2011
%                 
%                 try
%                     s= readannot_mat(dirann, aname,'tm_qrs');%OCT2012
%                     s.tm_qrs.pos = round(s.tm_qrs.pos*heasig.freq/1000);
%                     ext_anot=s.tm_qrs.pos(s.tm_qrs.pos >= t(1)&s.tm_qrs.pos <= t(2));%OCT2012
%                 catch
%                     
%                     s= readannot_mat(dirann, aname,'qrs');%DEC2011
%                     s.qrs.pos = round(s.qrs.pos*heasig.freq/1000);
%                     ext_anot=s.qrs.pos(s.qrs.pos >= t(1)&s.qrs.pos <= t(2));%Jul2011 % Jbolea 01/12/11
%                 end
%                 %                 aux_t=t/heasig.freq/60/60; %%%???? ver heasig.freq % Juan
%                 %                 28/03/2012
%                 
%                 %s=load ([dirann aname]);
%                 %kk=getfield(s);%JUL2011
%                 % eval(['ext_anot=s.qrs;']);
%             else error('QRS annotation file not found');
%             end
%         otherwise
%             messages.warnings=[messages.warnings {'Bad annotation format for external QRS annotation file; internal annotator would be considered instead.'}];
%             warning(char(messages.warnings(end)));
%             qrs_flag=0;
%     end
% end
% 
messages.setup.wavedet.qrs_flag=qrs_flag;
messages.setup.wavedet.leadsynth_flag=leadsynth_flag;
messages.setup.wavedet.flagaux=flagaux;

%%%%%%%%%%%%%%%%%%%%%%% Initiliazation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
inisamp = t(1);
endsamp = t(1)-1; %15MAR09
timeqrs1 = [];timeqrs2 = [];timeqrs3 = [];
lastqrs1 = 1;lastqrs2 = 1;lastqrs3 = 1;

if qrs_flag==0 || qrs_flag==2 % estimated maximum number of beats
    maxlength = round((t(2)-t(1))/messages.setup.wavedet.freq*3);
else maxlength=length(ext_anot);
end

if qrs_flag==1,
%     sel=find(ext_anot<messages.setup.wavedet.freq);
%     if ~isempty(sel)
%         ext_anot(sel)=[];
        maxlength= length(ext_anot);
%     end
end

nanvec = nan(1,maxlength);
position=struct('Pon',nanvec,'P',nanvec,'Poff',nanvec,...
    'Pprima',nanvec,'Pscale',nanvec,'Ptipo',nanvec,...
    'QRSon',nanvec,'Q',nanvec,'R',nanvec,'Rprima',nanvec,'S',nanvec,'QRSoff',nanvec,...
    'qrs',zeros(1,maxlength),...
    'R_inQRSoff',nanvec,'R_inQRSon',nanvec,...
    'Ton',nanvec,'T',nanvec,'Tprima', nanvec,'Toff',nanvec,'Ttipo',nanvec,...
    'Tscale',nanvec,'Ttipoon',nanvec,'Ttipooff',nanvec,'contadorToff',nanvec,...
    'QRSonsetcriteria',nanvec,'QRSoffcriteria',nanvec,...
    'QRSpa',nanvec,'QRSpp',nanvec,...
    'QRSmainpos',nanvec,'QRSmaininv',nanvec);


position1=position;
position2=position;
position3=position;
position0=position;
%messages.errors=[messages.errors errores.errores];


if qrs_flag==1, position.qrs=ext_anot;  end
% numlatdet = 0;      % Number of detected beats until now  % Juan
% 28/03/2012  no se usa despues
numlatdet1 = 0;
numlatdet2 = 0;
numlatdet3 = 0;

if ~isfield(messages.setup.wavedet,'filter_bank_design') || messages.setup.wavedet.filter_bank_design~=1 % Rute 24/11/11
    messages.setup.wavedet.filter_bank_design=0; % Rute 24/11/11
    [q1,q2,q3,q4,q5,messages] = qspfilt5(messages.setup.wavedet.freq,messages);  % WT equivalent filters % Rute 22/05/02
    if messages.status==1
        l1=length(q1);l2=length(q2);l3=length(q3);l4=length(q4);
        d1=floor((l1-1)/2);d2=floor((l2-1)/2);
        d3=floor((l3-1)/2);d4=floor((l4-1)/2);
        l5=length(q5);d5=floor((l5-1)/2); % Rute 22/05/02
        begoverlap = l5+2* messages.setup.wavedet.freq;    % l5 + 2 sec.
        endoverlap = d5;
%         messages.warnings=[messages.warnings {'Default filter banks are used.'}];
%         warning(char(messages.warnings(end)));
    else
        messages.setup.wavedet.filter_bank_design=1; % Rute 24/11/11
        messages.status=1;
    end
end
%%%% % Rute 24/11/11
if  messages.setup.wavedet.filter_bank_design==1;
    MaxScales=5;
    filters_cache_filename = ['wt_filters_' num2str(MaxScales) ' scales_' num2str(messages.setup.wavedet.freq) ' Hz.mat' ];
    % check it form setup
    if( exist(filters_cache_filename, 'file') )
        
        db_stat = dbstatus();
        bCaughtErrors = false;
        for ii = 1:length(db_stat)
            if( strcmpi(db_stat(ii).cond, 'caught error') )
                bCaughtErrors = true;
                dbclear if caught error
                break
            end
        end
        
        load( filters_cache_filename );
        
        if(bCaughtErrors)
            dbstop if caught error
        end
    else
        q_filters = qs_filter_design(1:MaxScales, messages.setup.wavedet.freq); %ver spf
        try
            save([fileparts(which('wavedet_3D')) sep filters_cache_filename],'q_filters')
        catch me
            me.messages = ['Fatal error saving ' filters_cache_filename ' on ' fileparts(which('wavedet_3D')) '; please check folder permissions'];
            messages.warnings = [messages.warnings {me.messages}];
        end
    end
    begoverlap=messages.setup.wavedet.freq;  % Juan, Julia 15/06/2012
    endoverlap=messages.setup.wavedet.freq;  % Juan, Julia 15/06/2012
    messages.warnings=[messages.warnings {'arbitrary sampling frequency designed filter banks are used.'}];
%     warning(char(messages.warnings(end)));
end
%%%%%%%%%%

% The two seconds overlap makes sure that the last/first beat
% are completely within the excerpt processed

samp=1; %heasig.spf_ecg*heasig.freq;    %Initialization of samp (1 sec)
intervaloall=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
heasigreading=heasig; %11JUL2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%% excerpt processing %%%%%%%%%%%%%%%%%%%%%%%%%%%
while ((endsamp+1) < t(2))
    fileevolution = [(endsamp+1) t(2)];
    
    firstnewsamp = samp(end);   % last sample analyzed
    endsamp = min(inisamp + nsamp -1,t(2));
    yy=1; %#ok<NASGU>
    %%%%% Protection 25/08/2010
    inisamp = round(inisamp);
    endsamp = round(endsamp);
%     display(['File Evolution (analyzing): ' timestr(round(inisamp/messages.setup.wavedet.freq*1000))...
%         ' to ' timestr(round(endsamp/messages.setup.wavedet.freq*1000))]);
    nsamp = round(nsamp);
    
    %%%%%%%%%%%%%%%%%%%%%%% loading ECG signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isnan(ft)||ft==0||ft==1||ft==2||ft==4||ft>4
%         if ft==0||ft==4 || ft==1 %MIT , Lund, Mortara %21MAR2011
%             [sig,messages1,heasig] = readsignal(sigdir,ecgnr,[inisamp endsamp],ft,heasigreading,leadreading,messages); % read leadreading
%             if messages1.status==0
%                 messages.errors=[messages.errors messages1.errors];
%                 warning(char(messages.errors(end)))
%                 messages.errors_desc=[messages.errors_desc messages1.errors_desc];
%                 messages.status=0;
%                 return;
%             else
%                 messages.warnings=[messages.warnings messages1.warnings ];%Jul2011
%                 if ~isfield(messages.setup,'readsignal') || ~isfield(messages.setup.readsignal,'file')
%                     messages.setup.readsignal.file=messages1.setup.file;
%                 end
%                 if ~isfield(messages.setup.readsignal,'t')
%                     messages.setup.readsignal.t=messages1.setup.t;
%                 else
%                     messages.setup.readsignal.t(2)=messages1.setup.t(2);
%                 end
%                 if ~isfield(messages.setup.readsignal,'nsig'), messages.setup.readsignal.nsig=messages1.setup.nsig; end
%             end
%         else
%             sig=sig_all(inisamp-t(1)+1:min(endsamp-t(1)+1,size(sig_all,1)),:);
%         end
        
        sig=sig_all(inisamp-t(1)+1:min(endsamp-t(1)+1,size(sig_all,1)));
        
        
%         if size(sig,1)==length(leadreading) % Jbolea 21/03/2011
%             sig=sig';
%         end
% 
%         n_auxleads=find(auxleads>0); %#ok<NASGU>

%         if isfield(heasig,'units')
%             for auxlead=1:length(leadreading) %
%                 if strcmp(deblank(upper(heasig.units(leadreading(auxlead),:))),'MV')||strcmp(deblank(upper(heasig.units(leadreading(auxlead),:))),'MILLIVOLTS')||...
%                         strcmp(deblank(upper(heasig.units(leadreading(auxlead),:))),'MILIVOLTIOS')
%                     sig(leadreading(auxlead),:) = sig(leadreading(auxlead),:)*1e3; % conversion to microV
%                     messages.warnings=[messages.warnings {'signal units changed: converted from milivolts to microvolts'}];
%                     warning(char(messages.warnings(end)))
%                     heasig.units=(char(ones(heasig.nsig,1)*'microvolts')); % Jbolea 21/3/2011
%                 elseif strcmp(deblank(upper(heasig.units(leadreading(auxlead),:))),'nV')||strcmp(deblank(upper(heasig.units(leadreading(auxlead),:))),'NANOVOLTS')||...
%                     strcmp(deblank(upper(heasig.units(leadreading(auxlead),:))),'NANOVOLTIOS')
%                     sig(:,leadreading(auxlead)) = (sig(:,leadreading(auxlead)))/1e3; % conversion to microV
%                     messages.warnings=[messages.warnings {'signal units changed: converted from volts to microvolts'}];
%                     warning(char(messages.warnings(end)))
%                     heasig.units=(char(ones(heasig.nsig,1)*'microvolts'));
%                 elseif strcmp(deblank(upper(heasig.units(leadreading(auxlead),:))),'V')||strcmp(deblank(upper(heasig.units(leadreading(auxlead),:))),'VOLTS')||...
%                         strcmp(deblank(upper(heasig.units(leadreading(auxlead),:))),'VOLTIOS')
%                     sig(leadreading(auxlead),:) = (sig(leadreading(auxlead),:))*1e6; % conversion to microV
%                     messages.warnings=[messages.warnings {'signal units changed: converted from volts to microvolts'}];
%                     warning(char(messages.warnings(end)))
%                     heasig.units=(char(ones(heasig.nsig,1)*'microvolts'));
%                 end
%             end
%         else
%             for auxlead=1:length(leadreading) %
%                 sig(leadreading(auxlead),:) = sig(leadreading(auxlead),:).* 1e3;  %% by default
%             end
%             messages.warnings=[messages.warnings {'no signal units given: milivolts assumed'}];
%             warning(char(messages.warnings(end)))
%             messages.warnings=[messages.warnings {'signal units changed: converted from milivolts to microvolts'}];
%             warning(char(messages.warnings(end)))
%             heasig.units=(char(ones(heasig.nsig,1)*'microvolts'));
%         end
        
%         if leadsynth_flag==1||leadsynth_flag==2
%             
%             if isnan(auxleads(8)) %lead II is not available
%                 sinal= leadcalc([sig(:,auxleads(1:7))  sig(:,auxleads(7))+sig(:,auxleads(9)) ],leadstring );  %Juan 23/06/2011
%             else
%                 sinal= leadcalc(sig,leadstring);
%             end
%             sig=sinal(:,~isnan(sinal(1,:)))';
%         elseif leadsynth_flag==3||leadsynth_flag==4||leadsynth_flag==5 % PC based delineation
%             if leadsynth_flag==4
%                 PCleads=1:size(sig,2); % considering all available leads
%             elseif leadsynth_flag==3 % considering 8 uncorrelated leads
%                 if ~isnan(auxleads(1:7)) & (~isnan(auxleads(8))||~isnan(auxleads(9))) %#ok<AND2>
%                     if isnan(auxleads(8)) %lead II is not available
%                         PCleads=[1:7 9]; %leads defined in leadreading, plus constructed lead II
%                         sig(:,end+1)=sig(:,7)+sig(:,8); %#ok<AGROW>
%                     else
%                         PCleads=1:8; %leads are already defined in leadreading
%                     end
%                 else
%                     messages.warnings=[messages.warnings {'There are not 8 uncorrelated leads, all available leads considered instead'}];
%                     warning(char(position.warningss(end)));
%                     PCleads=1:size(sig,2);
%                     leadsynth_flag=4;
%                 end
%             elseif leadsynth_flag==5
%                 if length(lead)>1 % considering precordial uncorrelated leads
%                     PCleads=2:7; %leads defined in leadreading, plus aVF
%                 else
%                     PCleads=[];
%                 end
%             end
%             if ~exist('PC_aux_ind','var') || isempty(PC_aux_ind)
%                 PC_aux_indaux= 1:length(sig);
%             else
%                 if flagaux == 1                    
%                     PC_aux_indaux=PC_aux_ind(PC_aux_ind>=firstnewsamp & PC_aux_ind<=endsamp)-firstnewsamp+1;
%                 end
%             end
%             %no NaN allowed in sig 29JUN2011
%             indx = find(isnan(sig(:,1)) == 1);
%             if ~isempty(indx)
%                 sig(indx,:) = sig(indx-1,:);
%             end
%             
%             [eigvectors,eigvalues]=eig(cov(sig(PC_aux_indaux,PCleads)));
%             [s_eigvalues,order]=sort(diag(eigvalues)); %#ok<ASGLU>
%             %pc=(eigvectors(:, order(end:-1:(end-length(lead)+1))))'*sig(:,PCleads)';
%             pc=(eigvectors(:, order(end:-1:(end-2))))'*sig(:,PCleads)'; %11JUL2011
%             if leadsynth_flag==5
%                 if isnan(auxleads(12)) %lead aVF is not available
%                     if ~isnan(auxleads(7)) && ~isnan(auxleads(8))
%                         sig=[sig(:,auxleads(8))-1/2*sig(:,auxleads(7)) pc(1:min(2,(length(lead)-1)),:)'];
%                     elseif ~isnan(auxleads(7)) && ~isnan(auxleads(9))
%                         sig=[sig(:,auxleads(9))+1/2*sig(:,auxleads(7)) pc(1:min(2,(length(lead)-1)),:)'];
%                     elseif ~isnan(auxleads(8)) && ~isnan(auxleads(9))
%                         sig=[1/2*(sig(:,auxleads(8))+sig(:,auxleads(9))) pc(1:min(2,(length(lead)-1)),:)'];
%                     else
%                         lead=auxleads(find(~isnan(auxleads([7:end 1:6])))); %#ok<FNDSB>
%                         messages.errors=[messages.errors {'lead aVF not available, PC based on precordial leads used only'}];
%                         warning(char(messages.errors(end)));
%                         sig= pc(1:min(3,(length(lead)-1)),:)';
%                         leadsynth_flag=4;
%                     end
%                 elseif length(lead)>1
%                     sig=[sig(:,7) pc'];
%                     sig=sig(:,lead);
%                 end
%             else
%                 sig=pc(lead,:)';
%             end
%         else %leadsynth_flag==0: should be leadsreading=lead
%             if lead<=size(sig,2)
%                 sig=sig(:,lead); %leadinfile=leadinfile(lead);
%                 try
%                     messages.warnings=[messages.warnings {[ 'the lead(s) chosen to be processed are the lead(s)' num2str(leadreading(lead))]}];
%                 catch me
%                     me.message = 'Error including warning information';
%                 end
%             else % if leads are in leadreading
%                 [C,IA,IB]=intersect(leadreading,lead); %#ok<NASGU>
%                 if length(C)==length(lead)
%                     sig=sig(:,IA);
%                     messages.warnings=[messages.warnings {['the lead(s) chosen to be processed are the lead(s) ' num2str(leadreading(IA))]}];
%                     lead=leadreading(IA); %7ABRIL2011
%                 else  %should not append
%                     lead=leadreading(1:min(3,end));
%                     messages.warnings=[messages.warnings {[ 'unable to interpret the lead(s) number to process (' num2str(lead) '): lead(s) ' num2str(leadreading) ' considered instead.' ]}];
%                 end
%             end
%         end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Modificado Juan 9/03/11 Rute 9May2011
    if size(sig,2)>1
        str = '';
    else
        str = '1';
    end
    %%%%%%%%%%%%%%%%%%%%%%% WT construction and parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if  messages.setup.wavedet.filter_bank_design==0;% RUTE 24NOV2011
        wt = wavt5(sig(:,1)',q1,q2,q3,q4,q5);  % WT equivalent filtering % 22/05/02 Rute including scale 5
        wt=wt(l5:end,1:5);            % Remove "incorrect" samples
        samp = (inisamp + l5 -1):endsamp-d5;
        % First l5-1 samples are not correctly filtered (border effect)
        % Last d5 samples are discarded in order to alineate all the
        % filtered signals, taking into acount the filter delays
        
        % Synchronizing filtered signals at different scales
        w1 = zeros(size(wt,1)-d5,5);
        w1(:,5) = wt(d5+1:end,5);
        w1(:,4) = wt(d4+1:end+d4-d5,4);
        w1(:,3) = wt(d3+1:end+d3-d5,3);
        w1(:,2) = wt(d2+1:end+d2-d5,2);
        w1(:,1) = wt(d1+1:end+d1-d5,1);
        
    else % RUTE 24NOV2011
        wtECG = qs_wt(sig,1:5, messages.setup.wavedet.freq, q_filters);% RUTE 24NOV2011 ver pf
        clear w1
        w1(1:size(wtECG,1),1:5) = wtECG(:,1,:);% RUTE 24NOV2011
        samp =inisamp:endsamp;% RUTE 24NOV2011
        l5=1;d5=0;
    end
    interval_to_processed=[inisamp+l5+1 endsamp-d5];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%check quality of the signal
    swt=w1.^2;
    % Initialization of the thresholds for QRS detection (eps)
    %step= messages.setup.wavedet.freq; %internal
%     if messages.setup.wavedet.sig_quality==2 %internal
%         step=heasig.spf_ecg*heasig.freq; %internal
%         indexes = signaltest(w1(:,1),step,1,sqrt(mean(swt(:,1))));
%         indexes(indexes>length(w1))=[];
%         
%         if ~isempty(indexes)
%             segments= indexes([find(diff(indexes)~=1) 1+find(diff(indexes)~=1)])';
%             segments=[indexes(1) segments(:)' indexes(end)];
%             segments=reshape(segments,2,length(segments)/2)'-d5-l5+1+t(1)+fileevolution(1)-2;
%             nsegments=size(segments,1);
%             resultados2xls(char(ones(nsegments,1)*ecgnr),[ones(nsegments,1)*[ft lead(1)] segments],[dirann_sig_quality aname_sig_quality '_sig_quality'],date);
%         end
%     elseif messages.setup.wavedet.sig_quality==1 %external
%         if ~iscell(aname_sig_quality)
%             aname_sig_quality_cell{1}=aname_sig_quality(:,1);
%         end
%         anots = readannot_mat(dirann_sig_quality,aname_sig_quality_cell{1},'no_quality_segments');
%         data = [anots.no_quality_segments.pos anots.no_quality_segments.info{1}];% timing and duration? units?
%         if  strcmpi(anots.no_quality_segments.units,'sec') ||   strcmpi(anots.no_quality_segments.units,'seconds')
%             data=messages.setup.wavedet.freq*data;
%         elseif  strcmpi(anots.no_quality_segments.units,'msec') || strcmpi(anots.no_quality_segments.units,'miliseconds')
%             data=messages.setup.wavedet.freq*data/1000;
%         else
%             messages.warnings=[messages.warnings {'no quality segments assumed to be in samples'}];
%         end
%         % for data in samples
%         indexes=data(data(:,1)>=interval_to_processed(1) & data(:,2)<=interval_to_processed(2),:);
%         indexes = intervalmatrix(indexes);
%         indexes(indexes>length(w1))=[];
%     else %no quality check
%         indexes=[];
%     end
%     if ~isempty(indexes)
%         w1(indexes,1)=0;
%         swt(indexes,1)=0;
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     eps= 0.5*sqrt(nanmean(swt));
    eps= rowvec(0.5*sqrt(nanmedian(swt))) .* rowvec(messages.setup.wavedet.QRS_detection_thr);
    eps(4) = eps(4)*2;
    
    %%%%%%% multilead%%% Rute 02.Dec.04
%     if size(sig,2)>1
%         % lead 2
%         if  messages.setup.wavedet.filter_bank_design==0;% RUTE 24NOV2011
%             wt = wavt5(sig(:,2)',q1,q2,q3,q4,q5);  % WT equivalent filtering % 22/05/02 Rute including scale 5
%             wt=wt(l5:end,1:5);            % Remove "incorrect" samples
%             %samp = (inisamp + l5 -1):endsamp-d5;
%             % First l5-1 samples are not correctly filtered (border effect)
%             % Last d5 samples are discarded in order to alineate all the
%             % filtered signals, taking into acount the filter delays
%             % Synchronizing filtered signals at different scales
%             w2 = zeros(size(wt,1)-d5,5);
%             w2(:,5) = wt(d5+1:end,5);
%             w2(:,4) = wt(d4+1:end+d4-d5,4);
%             w2(:,3) = wt(d3+1:end+d3-d5,3);
%             w2(:,2) = wt(d2+1:end+d2-d5,2);
%             w2(:,1) = wt(d1+1:end+d1-d5,1);
%         else
%             clear w2
%             w2(1:size(wtECG,1),1:5) = wtECG(:,2,:);% RUTE 24NOV2011
%         end
%         swt=w2.^2;
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%check quality of the signal
%         if messages.setup.wavedet.sig_quality==1 %external
%             if iscell(aname_sig_quality) && size(aname_sig_quality,2)>1
%                 aname_sig_quality_cell{2}=aname_sig_quality(:,2);
%                 anots = readannot_mat(dirann_sig_quality,aname_sig_quality_cell{2},'no_quality_segments');
%                 data = [anots.no_quality_segments.pos anots.no_quality_segments.info{1}];% timing and duration? units?
%                 if  strcmpi(anots.no_quality_segments.units,'sec') ||   strcmpi(anots.no_quality_segments.units,'seconds')
%                     data=messages.setup.wavedet.freq*data;
%                 elseif  strcmpi(anots.no_quality_segments.units,'msec') || strcmpi(anots.no_quality_segments.units,'miliseconds')
%                     data=messages.setup.wavedet.freq*data/1000;
%                 else
%                     messages.warnings=[messages.warnings {'no quality segments assumed to be in samples'}];
%                 end                % for data in samples
%                 indexes=data(data(:,1)>=interval_to_processed(1) & data(:,2)<=interval_to_processed(2),:);
%                 indexes = intervalmatrix(indexes);
%             end
%             indexes(indexes>length(w2))=[];
%         else %no quality check
%             indexes=[];
%         end
%         if ~isempty(indexes)
%             w2(indexes,1)=0;
%             swt(indexes,1)=0;
%         end
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         % Initialization of the thresholds for QRS detection (eps)
%         eps(2,:)= 0.5*sqrt(mean(swt));
%         eps(2,4) = eps(2,4)*2;
%         
%         if size(sig,2)>2
%             % lead 3
%             if  messages.setup.wavedet.filter_bank_design==0;% RUTE 24NOV2011
%                 wt = wavt5(sig(:,3)',q1,q2,q3,q4,q5);  % WT equivalent filtering % 22/05/02 Rute including scale 5
%                 wt=wt(l5:end,1:5);            % Remove "incorrect" samples
%                 %samp = (inisamp + l5 -1):endsamp-d5;
%                 % First l5-1 samples are not correctly filtered (border effect)
%                 % Last d5 samples are discarded in order to alineate all the
%                 % filtered signals, taking into acount the filter delays
%                 
%                 % Synchronizing filtered signals at different scales
%                 w3 = zeros(size(wt,1)-d5,5);
%                 w3(:,5) = wt(d5+1:end,5);
%                 w3(:,4) = wt(d4+1:end+d4-d5,4);
%                 w3(:,3) = wt(d3+1:end+d3-d5,3);
%                 w3(:,2) = wt(d2+1:end+d2-d5,2);
%                 w3(:,1) = wt(d1+1:end+d1-d5,1);
%             else
%                 clear w3
%                 w3(1:size(wtECG,1),1:5) = wtECG(:,3,:);% RUTE 24NOV2011
%             end
%             swt=w3.^2;
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%check quality of the signal
%             if messages.setup.wavedet.sig_quality==1 %external
%                 if iscell(aname_sig_quality) && size(aname_sig_quality,2)>2
%                     aname_sig_quality_cell{3}=aname_sig_quality(:,3);
%                     anots = readannot_mat(dirann_sig_quality,aname_sig_quality_cell{3},'no_quality_segments');
%                     data = [anots.no_quality_segments.pos anots.no_quality_segments.info{1}];% timing and duration? units?
%                     if  strcmpi(anots.no_quality_segments.units,'sec') ||   strcmpi(anots.no_quality_segments.units,'seconds')
%                         data=messages.setup.wavedet.freq*data;
%                     elseif  strcmpi(anots.no_quality_segments.units,'msec') || strcmpi(anots.no_quality_segments.units,'miliseconds')
%                         data=messages.setup.wavedet.freq*data/1000;
%                     else
%                         messages.warnings=[messages.warnings {'no quality segments assumed to be in samples'}];
%                     end
%                     % for data in samples
%                     indexes=data(data(:,1)>=interval_to_processed(1) & data(:,2)<=interval_to_processed(2),:);
%                     indexes = intervalmatrix(indexes);
%                 end
%                 indexes(indexes>length(w3))=[];
%             else %no quality check
%                 indexes=[];
%             end
%             if ~isempty(indexes)
%                 w3(indexes,1)=0;
%                 swt(indexes,1)=0;
%             end
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % Initialization of the thresholds for QRS detection (eps)
%             eps(3,:)= 0.5*sqrt(mean(swt));
%             eps(3,4) = eps(3,4)*2;
%         end
%     end
%     
    sig=sig(l5:end-d5,:); % piece of signal processed in this iteration
    %messages.lixo=sig;
    eps=eps';
    clear wt swt;
    
%%    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% beat detection %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if qrs_flag==0 || qrs_flag==2   % The QRS fiducial point is calculated with wavedet
        messages.sig=sig;
        [position1,timeqrs1,lastqrs1,intervalo1,numlatdet1,time1,messages]= fiducialf(position1,firstnewsamp,samp,heasig,w1,eps,1,lastqrs1,timeqrs1,numlatdet1,ultimo_anot,messages);
        timeqrs=timeqrs1; % comentar isto pode dar erro no single lead
        indexes=1:length(timeqrs1);
        
%         if size(sig,2)>1
%             [position2,timeqrs2,lastqrs2,intervalo2,numlatdet2,time2,messages]= fiducialf(position2,firstnewsamp,samp,heasig,w2,eps,2,lastqrs2,timeqrs2,numlatdet2,ultimo_anot,messages);
%             if size(sig,2)>2
%                 [position3,timeqrs3,lastqrs3,intervalo3,numlatdet3,time3,messages]= fiducialf(position3,firstnewsamp,samp,heasig,w3,eps,3,lastqrs3,timeqrs3,numlatdet3,ultimo_anot,messages);
%             else
%                 timeqrs3=[];
%             end
%             [timeqrs,indexes]=qrscandidatesnew(timeqrs1,timeqrs2,timeqrs3,heasig,messages.setup.wavedet.refrper(end)/2* messages.setup.wavedet.freq);
%             
%             %2JUN08 %RUTE no candidates at all do not mark
%             if ~isempty(timeqrs)
%                 lastqrs1=timeqrs(1,~isnan(timeqrs(1,:)));
%                 if ~isempty(lastqrs1)
%                     lastqrs1=lastqrs1(end);
%                 end
%                 
%                 if sum(~isnan(timeqrs(2,:)))>0
%                     lastqrs2=timeqrs(2,~isnan(timeqrs(2,:)));
%                 end
%                 if ~isempty(lastqrs2)
%                     lastqrs2=lastqrs2(end);
%                 end
%                 if size(sig,2)>2
%                     if sum(~isnan(timeqrs(3,:)))>0
%                         lastqrs3=timeqrs(3,~isnan(timeqrs(3,:)));
%                     end
%                     if ~isempty(lastqrs3)
%                         lastqrs3=lastqrs3(end);
%                     end
%                 end
%             end
%         end
        
        intervaloall=[intervaloall(end)+1 size(timeqrs,2)];
        intervalo=intervaloall;
        time= timeqrs(:,intervaloall(1):intervaloall(end))-samp(1)+1;
        
        if ~isempty(time) && size(sig,2)>1 %16OUT08
            position0.qrs(intervalo(1):intervalo(end))=round(nanmedian(time)+samp(1)-1); % 09AGO2011
            position.qrs(intervalo(1):intervalo(end))=round(nanmedian(time)+samp(1)-1);  % 09AGO2011
        end
    else           % The QRS fiducial point is read from external annotator % multilead need to be checked
        
        
        first = max(firstnewsamp-ceil(messages.setup.wavedet.freq), samp(1)-1+ceil( messages.setup.wavedet.freq*0.050));
        % first is calculated according with the first lines of fiducial JP
        sel=find(ext_anot>=first & ext_anot<=samp(end)); %
%         position1.qrs(sel)=ext_anot(sel); %4ABRIL2011
        position1.qrs = rowvec(ext_anot);
        timeqrs(1,sel)=position1.qrs(sel)-samp(1)+1; %global
        time=timeqrs(1,sel);%4ABRIL2011 ?????%local
        indexes= 1:length(timeqrs(1,:));%4ABRIL2011
        lastqrs1=timeqrs(1,~isnan(timeqrs(1,:)));
        if ~isempty(lastqrs1)
            lastqrs1=lastqrs1(end);
        end
        if ~isempty(sel)
            intervaloall=[sel(1) sel(end)];
            % time=[timeqrs(:,intervaloall(1):intervaloall(end))-samp(1)+1];
        else
            intervaloall=[];
            %time=[];
        end
        intervalo=intervaloall;
        intervalo1=intervalo;%4ABRIL2011
        time1=time;%4ABRIL2011
        
        if size(sig,2)>1
            position2.qrs(sel)=ext_anot(sel);%4ABRIL2011
            timeqrs(2,sel)=position2.qrs(sel)-samp(1)+1;
            time(2,:)= time(1,:);%4ABRIL2011
            indexes=[indexes; indexes];%#ok<AGROW>%4ABRIL2011
            if sum(~isnan(timeqrs(2,:)))>0
                lastqrs2=timeqrs(2,~isnan(timeqrs(2,:)));
            end
            if ~isempty(lastqrs2)
                lastqrs2=lastqrs2(end);
            end
            intervalo2=intervalo;%4ABRIL2011
            time2=time;
            
            if size(sig,2)>2
                position3.qrs(sel)=ext_anot(sel);%4ABRIL2011
                timeqrs(3,sel)=position3.qrs(sel)-samp(1)+1;
                time=[time;time(1,:)];%#ok<AGROW> %4ABRIL2011
                indexes=[indexes; indexes];%#ok<AGROW>%4ABRIL2011
                if sum(~isnan(timeqrs(3,:)))>0
                    lastqrs3=timeqrs(3,~isnan(timeqrs(3,:)));
                end
                if ~isempty(lastqrs3)
                    lastqrs3=lastqrs3(end);
                end
                
                intervalo3=intervalo;%4ABRIL2011
                time3=time;
            else
                indexes=[indexes;NaN*(1:length(timeqrs(1,:)))];%#ok<AGROW>%9MAYO2011
            end
        end
        
        if ~isempty(time) &&  size(sig,2)>2 %Rute 13Set06 %by default position is the median mark
            position0.qrs(intervalo(1):intervalo(end))=round(nanmedian(time)+samp(1)-1); % 09AGO2011
            position.qrs(intervalo(1):intervalo(end))=round(nanmedian(time)+samp(1)-1); % 09AGO2011
        end
        
        % If the P-wave is not in the present segment, discard the beat, because the
        % the beat should have been detected in the last segment.
        
        %%%%%%%%%%%%%%%%% Rute 03.Dec.04?
        %         if abs(samp(time(1)) - lastqrs)< 0.275*heasig.spf_ecg*heasig.freq,
        %             time(1)=[];  sel(1)=[];              % In order not to repeat beats
        %         end
        %         if (length(sig)-time(end) < 0.45*heasig.spf_ecg*heasig.freq),
        %             time(end)=[]; sel(end)=[];
        %         end
        % If last beat's T-wave is not in 'sig', the beat is detected in next one
        %%% timeqrs = [timeqrs samp(time)];  % QRS times ??????????????? Rute 03.Dec.04?
        %lastqrs = samp(time(end));Rute 03.Dec.04?
        % numeration of the beats processed in this segment
    end
    if size(sig,2)>1
        timenew=round(nanmedian(time)); %RUTE 09/07/2007
    else
        timenew=time;
    end
    
    
    if( ~messages.setup.wavedet.QRS_detection_only )
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% QRS delineation %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %if ~isempty(time(1,:)) %2JUN08 RUTE
        if ~isempty(time) %2JUN08 RUTE
            QRSon=indexes(:,intervalo(1):intervalo(end));
            qrspiconsetall=NaN*ones(3,intervalo(end)-intervalo(1)+1);
            QRSoff=indexes(:,intervalo(1):intervalo(end));
            qrspicoffsetall=NaN*ones(3,intervalo(end)-intervalo(1)+1);
            if ~isempty(intervalo1)
                [position1,qrspiconset1,qrspicoffset1,messages]= qrswavef(heasig,samp,time1,position1,w1,intervalo1,messages);
                auxindexes=indexes(1,intervalo(1):intervalo(end));
                %NOTE: CHANGED 10.SET.05
                qrspiconsetall(1,(~isnan(indexes(1,intervalo(1):intervalo(end)))&(QRSon(1,:)>=intervalo1(1)&QRSon(1,:)<=intervalo1(2))))=qrspiconset1(QRSon(1,~isnan(QRSon(1,:))&(QRSon(1,:)>=intervalo1(1)&QRSon(1,:)<=intervalo1(2)))-intervalo1(1)+1); %17Mar06
                QRSon(1,~isnan(QRSon(1,:)))=position1.QRSon(QRSon(1,~isnan(QRSon(1,:))));

                qrspicoffsetall(1,(~isnan(indexes(1,intervalo(1):intervalo(end)))&(QRSoff(1,:)>=intervalo1(1)&QRSoff(1,:)<=intervalo1(2))))=qrspicoffset1(QRSoff(1,~isnan(QRSoff(1,:))&(QRSoff(1,:)>=intervalo1(1)&QRSoff(1,:)<=intervalo1(2)))-intervalo1(1)+1);
                QRSoff(1,~isnan(QRSoff(1,:)))=position1.QRSoff(QRSoff(1,~isnan(QRSoff(1,:))));
            end
            if size(sig,2)>1
                [position2,qrspiconset2,qrspicoffset2,messages]= qrswavef(heasig,samp,time2,position2,w2,intervalo2,messages);
                if ~isempty(intervalo2)
                    auxindexes=indexes(2,intervalo(1):intervalo(end));
                    qrspiconsetall(2,(~isnan(indexes(2,intervalo(1):intervalo(end)))&(QRSon(2,:)<=intervalo2(2)&QRSon(2,:)>=intervalo2(1))))=qrspiconset2(QRSon(2,~isnan(QRSon(2,:))&(QRSon(2,:)<=intervalo2(2)&QRSon(2,:)>=intervalo2(1)))-intervalo2(1)+1); %Rute 17Mar06
                    QRSon(2,~isnan(QRSon(2,:)))=position2.QRSon(QRSon(2,~isnan(QRSon(2,:))));
                    qrspicoffsetall(2,(~isnan(indexes(2,intervalo(1):intervalo(end)))&(QRSoff(2,:)<=intervalo2(2)&QRSoff(2,:)>=intervalo2(1))))=qrspicoffset2(QRSoff(2,~isnan(QRSoff(2,:))&(QRSoff(2,:)<=intervalo2(2)&QRSoff(2,:)>=intervalo2(1)))-intervalo2(1)+1); %Rute 17Mar06
                    QRSoff(2,~isnan(QRSoff(2,:)))=position2.QRSoff(QRSoff(2,~isnan(QRSoff(2,:))));
                end
                if size(sig,2)>2 && ~isempty(intervalo3)
                    [position3,qrspiconset3,qrspicoffset3,messages]= qrswavef(heasig,samp,time3,position3,w3,intervalo3,messages); %#ok<ASGLU>
                    auxindexes=indexes(3,intervalo(1):intervalo(end));
                    if max(QRSon(3,~isnan(QRSon(3,:))&(QRSon(3,:)>=intervalo3(1)&QRSon(3,:)<=intervalo3(2)))-intervalo3(1)+1)> length(qrspiconset3)%17.Mar.06
                        qrspiconsetall(3,(QRSon(3,~isnan(QRSon(3,:))&(QRSon(3,:)<=intervalo3(2)&QRSon(3,:)>=intervalo3(1)))-intervalo3(1)+1)<=length(qrspiconsetall(3,(~isnan(indexes(3,intervalo(1):intervalo(end)))&(QRSon(3,:)<=intervalo3(2)&QRSon(3,:)>=intervalo3(1)))))&(~isnan(indexes(3,intervalo(1):intervalo(end)))&(QRSon(3,:)>=intervalo3(1))))=qrspiconset3((QRSon(3,~isnan(QRSon(3,:))&(QRSon(3,:)>=intervalo3(1)))-intervalo3(1)+1)<=length(qrspiconsetall(3,(~isnan(indexes(3,intervalo(1):intervalo(end)))&(QRSon(3,:)>=intervalo3(1)))))& QRSon(3,~isnan(QRSon(3,:))&(QRSon(3,:)>=intervalo3(1)))-intervalo3(1)+1);%17Mar06
                    else%16.Fev.06
                        qrspiconsetall(3,(~isnan(indexes(3,intervalo(1):intervalo(end)))&(QRSon(3,:)>=intervalo3(1))))=qrspiconset3(QRSon(3,~isnan(QRSon(3,:))&(QRSon(3,:)>=intervalo3(1)))-intervalo3(1)+1);
                    end %16.Fev.06
                    QRSon(3,~isnan(QRSon(3,:)))=position3.QRSon(QRSon(3,~isnan(QRSon(3,:))));

                    if max(QRSoff(3,~isnan(QRSoff(3,:))&(QRSoff(3,:)>=intervalo3(1)&QRSoff(3,:)<=intervalo3(2)))-intervalo3(1)+1)> length(qrspiconset3)%17.Mar.06
                        qrspicoffsetall(3,(QRSoff(3,~isnan(QRSoff(3,:))&(QRSoff(3,:)<=intervalo3(2)&QRSoff(3,:)>=intervalo3(1)))-intervalo3(1)+1)<=length(qrspicoffsetall(3,(~isnan(indexes(3,intervalo(1):intervalo(end)))&(QRSoff(3,:)<=intervalo3(2)&QRSoff(3,:)>=intervalo3(1)))))&(~isnan(indexes(3,intervalo(1):intervalo(end)))&(QRSoff(3,:)>=intervalo3(1))))=qrspiconset3((QRSoff(3,~isnan(QRSoff(3,:))&(QRSoff(3,:)>=intervalo3(1)))-intervalo3(1)+1)<=length(qrspicoffsetall(3,(~isnan(indexes(3,intervalo(1):intervalo(end)))&(QRSoff(3,:)>=intervalo3(1)))))& QRSoff(3,~isnan(QRSoff(3,:))&(QRSoff(3,:)>=intervalo3(1)))-intervalo3(1)+1);%17Mar06
                    else%16.Fev.06
                        qrspicoffsetall(3,(~isnan(indexes(3,intervalo(1):intervalo(end)))&(QRSoff(3,:)>=intervalo3(1))))=qrspiconset3(QRSoff(3,~isnan(QRSoff(3,:))&(QRSoff(3,:)>=intervalo3(1)))-intervalo3(1)+1);
                    end %16.Fev.06
                    QRSoff(3,~isnan(QRSoff(3,:)))=position3.QRSoff(QRSoff(3,~isnan(QRSoff(3,:))));
                end
            end
            auxqrspiconset=nanmin(qrspiconsetall);
            auxqrspicoffset=nanmax(qrspicoffsetall);
            aux=find(isnan(auxqrspiconset));
            auxaux=intervalo(1):intervalo(end);
            %
            if ~isempty(aux)  % one beat is from the former interval  1.AGO.07
%                 if aux(1)>1
%                     messages.errors=[messages.errors {['QRS onset not found at beat(s): ' num2str(auxaux(aux ))]}];
%                     warning(char(messages.errors(end)));
%                 end
                auxqrspiconset(aux)=[];
                QRSon(:,aux)=[];  %31.07.07
                auxqrspicoffset(aux)=[];
                QRSoff(:,aux)=[]; %31.07.07
                intervalo(end)=intervalo(end)-aux(end); %16MAR09
                time(:,aux)=[];%31.07.07
                timenew(:,aux)=[];%31.07.07
                if length(auxaux)>size(QRSon,2)
                    auxaux(aux)=[];
                end
            end

            if size(sig,2)>1 %Rute 16Out08
                %by default QRS onset in the earliest of the marks
                %by default QRS end in the latest of the marks
                position0.QRSon(auxaux)=min(QRSon);
                position.QRSon(auxaux)=min(QRSon);
                position0.QRSoff(auxaux)=max(QRSoff);
                position.QRSoff(auxaux)=max(QRSoff);
            end
            if size(sig,2)>1
                if ~exist('w3','var')
                    w3=[];
                end
                timeaux=time; %time is sorted %18.07.07
                auxaux=NaN*ones(3,intervalo(end)-intervalo(1)+1);
                auxindexes=indexes(:,intervalo(1):intervalo(end));
                auxaux(1,~isnan(indexes(1,intervalo(1):intervalo(end))))=position1.R(auxindexes(1,~isnan(indexes(1,intervalo(1):intervalo(end)))));
                auxaux(2,~isnan(indexes(2,intervalo(1):intervalo(end))))=position2.R(auxindexes(2,~isnan(indexes(2,intervalo(1):intervalo(end)))));
                auxaux(3,~isnan(indexes(3,intervalo(1):intervalo(end))))=position3.R(auxindexes(3,~isnan(indexes(3,intervalo(1):intervalo(end)))));
                position.R(intervalo(1):intervalo(end))=nanmin(auxaux);

                %timenew is the median, use min instead?
                [position0,position,messages]=multiqrson(heasig,w1,w2,w3,auxqrspiconset,QRSon,timeaux,timenew,position0,position,intervalo,samp,messages);
                timeaux=time;

                position.R(intervalo(1):intervalo(end))=round(nanmedian(auxaux)); %RUTE 09AGO2011 %MEDIAN %by default position is the median mark of all R
                %timenew is the median, use max instead?
                [position0,position,messages]=multiqrsend(heasig,w1,w2,w3,auxqrspicoffset,QRSoff,timeaux,timenew,position0,position, intervalo,samp,messages);
                position.R(intervalo(1):intervalo(end))=round(nanmedian(auxaux));%RUTE 09AGO2011
                %            squaresignal=(sig(:,1)/max(sig(:,1))).^2+(sig(:,2)/max(sig(:,2))).^2+(sig(:,3)/max(sig(:,3))).^2;
                %             for ii=intervalo(1):intervalo(end)
                %             if ~isempty(position.QRSon(ii):position.QRSoff(ii))
                %                 [M,position.qrsnew(ii)]=max(squaresignal(position.QRSon(ii):position.QRSoff(ii)));
                %                 position.qrsnew(ii)=position.qrsnew(ii)+position.QRSon(ii)-1;
                %             end
                %             end
            end
            %         lastqrs1=timeqrs(1,end);
            %%%%%%%%%%%%%%%% caso em que rr tem dimensao inferior ao necessario
            %%%%%%%%%%%%%%%% protection with respect to ultimo_anot

            while eval(['position' str '.QRSon(intervalo(1))<ultimo_anot'])
                if ultimo_anottyp==10; % Toff
                    if eval(['(position' str '.qrs(intervalo(1))-position' str '.QRSon(intervalo(1))) < (position' str '.Toff(intervalo(1)-1)-position' str '.T(intervalo(1)-1))'])
                        eval(['position' str '.Toff(intervalo(1)-1)=NaN;'])
                        if eval(['isnan(position' str '.Tprima(intervalo(1)-1))'])
                            eval(['ultimo_anot=position' str '.T(intervalo(1)-1);'])
                            ultimo_anottyp=8;
                        else
                            eval(['ultimo_anot=position' str '.Tprima(intervalo(1)-1);'])
                            ultimo_anottyp=9;
                        end
                    else
                        eval(['position' str '.QRSon(intervalo(1))=NaN;'])
                    end
                elseif ultimo_anottyp==9 || ultimo_anottyp==8 || ultimo_anottyp==5; % Tprima | Tpeak | qrs
                    eval(['position' str '.QRSon(intervalo(1))=NaN;'])
                else  % ultimo_anottyp==6 % QRSoff
                    if eval(['(position' str '.qrs(intervalo(1))-position' str '.QRSon(intervalo(1))) < (position' str '.QRSoff(intervalo(1)-1)-position' str '.qrs(intervalo(1)-1))'])
                        eval(['position' str '.qrsoff(intervalo(1)-1)=NaN;'])
                        eval(['ultimo_anot=position' str '.qrs(intervalo(1)-1);'])
                        ultimo_anottyp=5;
                    else
                        eval(['position' str '.QRSon(intervalo(1))=NaN;'])
                    end
                end
            end
            % Modificado Juan 9/03/11
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if qrs_flag~=2
                %%%%%%%%%%%%%%%%%%%%%%%%%%%% T wave
                if ~isempty(time1)
                    [position1,janelas1,messages]= twavef(heasig,samp,time1,position1,w1,intervalo1,messages); %21AGO09
                end
                %end
                if size(sig,2)>1 && ~isempty(time2)
                    if length(time2(1,:))>1%23SET08
                        [position2,janelas2,messages]= twavef(heasig,samp,time2,position2,w2,intervalo2,messages);%21AGO09
                    end
                    if size(sig,2)>2  && ~isempty(intervalo3)
                        if length(time3(1,:))>1%23SET08
                            [position3,janelas3,messages]= twavef(heasig,samp,time3,position3,w3,intervalo3,messages);%21AGO09
                        end
                    end
                end
                if size(sig,2)>1
                    subintervalo=indexes(:,intervalo(1):intervalo(end));
                    %NOTE are excluded the ones from the previous seg
                    if ~isempty(intervalo1)
                        aux1=subintervalo(1,~isnan(subintervalo(1,:)))-intervalo1(1)+1;
                    end
                    if size(sig,2)>1 & ~isempty(intervalo2) %#ok<AND2>
                        aux2=subintervalo(2,~isnan(subintervalo(2,:)))-intervalo2(1)+1;
                    end
                    if size(sig,2)>2 & ~isempty(intervalo3) %#ok<AND2>
                        aux3=subintervalo(3,~isnan(subintervalo(3,:)))-intervalo3(1)+1;
                    end
                    %26.04.05
                    janelas=NaN*ones(size(subintervalo,2),9);
                    janelas(:,1)=(intervalo(1):intervalo(end))';
                    janelas(:,4)=(intervalo(1):intervalo(end))';
                    janelas(:,7)=(intervalo(1):intervalo(end))';
                    if ~isempty(intervalo1)
                        janelas(~isnan(subintervalo(1,:)-intervalo1(1)+1),2)=aux1';
                    end
                    if ~isempty(intervalo2)
                        janelas(~isnan(subintervalo(2,:)-intervalo2(1)+1),5)=aux2';
                    end
                    if ~isempty(intervalo1)
                        aux1=aux1(aux1>0);
                    end
                    aux2=aux2(aux2>0);
                    if size(sig,2)>2 & ~isempty(intervalo3) %#ok<AND2>
                        janelas(~isnan(subintervalo(3,:)-intervalo3(1)+1),8)=aux3';
                        aux3=aux3(aux3>0);
                    end
                    janelas(janelas<=0)=NaN;
                    if exist('janelas1','var')&& sum(~isnan(janelas(:,2)))~=0
                        janelas(~isnan(janelas(:,2)) & janelas(:,2)>0,2:3)=janelas1(aux1,2:3);
                    end
                    if exist('janelas2','var')& ~isempty(janelas(~isnan(janelas(:,5)) & janelas(:,5)>0 ,5:6)) %#ok<AND2>
                        if max(aux2) > length(janelas2)%17.Mar.06
                            janelas(~isnan(janelas(:,5)) & janelas(:,5)>0 ,5:6)=[janelas2((aux2) <= length(janelas2),2:3); NaN*ones(sum(aux2 > length(janelas2)),2)];%14.Mar.06
                        else%16.Fev.06
                            janelas(~isnan(janelas(:,5)) & janelas(:,5)>0 ,5:6)=janelas2(aux2,2:3);
                        end%16.Fev.06
                    end
                    if size(sig,2)>2 & ~isempty(intervalo3) %#ok<AND2>
                        %ver altera�ao por erro de indexes em irec5 de
                        %politec!!!!%16.Fev.06
                        if exist('janelas3','var')
                            if max(aux3) > length(janelas3)%16.Fev.06
                                %janelas(~isnan(janelas(:,8)) & janelas(:,8)>0 ,8:9)=[janelas3(aux3(1:(end-1)),2:3);NaN NaN];%16.Fev.06
                                janelas(~isnan(janelas(:,8)) & janelas(:,8)>0 ,8:9)=[janelas3((aux3) <= length(janelas3),2:3); NaN*ones(sum(aux3 > length(janelas3)),2)];%14.Mar.06
                            else%16.Fev.06
                                janelas(~isnan(janelas(:,8)) & janelas(:,8)>0 ,8:9)=janelas3(aux3,2:3);
                            end%16.Fev.06
                        end
                    end
                    %cases were begwin>endwin!!!
                    janelas((janelas(:,3))<=janelas(:,2)|isnan(janelas(:,3)),2)=NaN;
                    janelas(janelas(:,3)<=janelas(:,2)|isnan(janelas(:,2)),3)=NaN;
                    janelas(janelas(:,6)<=janelas(:,5)|isnan(janelas(:,6)),5)=NaN;
                    janelas(janelas(:,6)<=janelas(:,5)|isnan(janelas(:,5)),6)=NaN;
                    janelas(janelas(:,9)<=janelas(:,8)|isnan(janelas(:,9)),8)=NaN;
                    janelas(janelas(:,9)<=janelas(:,8)|isnan(janelas(:,8)),9)=NaN;
                    begwin= min(janelas(:,[2 5 8])'); %#ok<UDIM>
                    endwin= max(janelas(:,[3 6 9])'); %#ok<UDIM>
                    intreg=[begwin' endwin'];

                    count=count+1;
                    if ~exist('w3','var')
                        w3=[];
                    end
                    [position0,position,messages]= delineate3D(w1,w2,w3,intreg,timenew,position0,position,intervalo,position1,position2,position3,indexes,samp,[],messages); %   05MAY2011
                end
                %end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 24AGO2011
                if ~exist('str','var'), str=''; end
                eval(['position' str '.QRSon(position' str '.R<=(position' str '.QRSon+1))=NaN;']);
                eval(['position' str '.QRSoff(position' str '.QRSoff>=(position' str '.Ton-1))=NaN;']);
                eval(['position' str '.QRSoff(position' str '.QRSoff>=(position' str '.T-1))=NaN;']);
                eval(['position' str '.QRSoff(position' str '.QRSoff<=(position' str '.qrs+1))=NaN;']);
                eval(['position' str '.QRSoff(position' str '.R>=(position' str '.QRSoff-1))=NaN;']);
                eval(['position' str '.Toff(position' str '.Toff<=(position' str '.T+1))=NaN;']);
                eval(['position' str '.Ton(position' str '.Ton>=(position' str '.T-1))=NaN;']);
                eval(['ii=1:sum(position' str '.qrs>0);']);
                eval(['position' str '.QRSoff(position' str '.QRSoff(ii(1:end-1))>=(position' str '.R(ii(2:end))))=NaN;']);
                eval(['position' str '.QRSoff(position' str '.QRSoff(ii(1:end-1))>=(position' str '.qrs(ii(2:end))))=NaN;']);
                eval(['position' str '.Toff(position' str '.Toff(ii(1:end-1))>=position' str '.qrs(ii(2:end)))=NaN;']);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%% P wave %%%%%%%%%%%%%%%%%%%%%%%%%%%
                %11.04.05 P wave delineation inside the condition of existing some
                [position1,picon_all,picoff_all,messages]= pwavef(heasig,samp,time1,position1,w1,intervalo1,ultimo_anot,messages); %#ok<ASGLU>
                if size(sig,2)>1
                    [position2,picon_all,picoff_all,messages]= pwavef(heasig,samp,time2,position2,w2,intervalo2,ultimo_anot,messages); %#ok<ASGLU>
                    if size(sig,2)>2 & ~isempty(intervalo3) %#ok<AND2>
                        [position3,picon_all,picoff_all,messages]= pwavef(heasig,samp,time3,position3,w3,intervalo3,ultimo_anot,messages); %#ok<ASGLU>
                        auxaux(1,~isnan(indexes(1,intervalo(1):intervalo(end))))=position1.P(auxindexes(1,~isnan(indexes(1,intervalo(1):intervalo(end)))));
                        auxaux(2,~isnan(indexes(2,intervalo(1):intervalo(end))))=position2.P(auxindexes(2,~isnan(indexes(2,intervalo(1):intervalo(end)))));
                        auxaux(3,~isnan(indexes(3,intervalo(1):intervalo(end))))=position3.P(auxindexes(3,~isnan(indexes(3,intervalo(1):intervalo(end)))));
                        position0.P(intervalo(1):intervalo(end))=round(nanmedian(auxaux)); %RUTE 09AGO2011; %MEDIAN %by default position is the median mark of all R
                        [position0,position,messages]= delineateP3D(heasig,w1,w2,w3,timenew,position0,position,intervalo,samp,ultimo_anot,messages);
                    end
                end
                %check the position structure size
                eval(['position' str '.Pon((length(position' str '.Pon)+1):length(position' str '.qrs))=NaN;']); %6May2011
                eval(['position' str '.P((length(position' str '.P)+1):length(position' str '.qrs))=NaN;']); %6May2011
                eval(['position' str '.Poff((length(position' str '.Poff)+1):length(position' str '.qrs))=NaN;']); %6May2011
                eval(['position' str '.Poff((length(position' str '.Poff)+1):length(position' str '.qrs))=NaN;']); %6May2011
                if size(sig,2)>2 %6May2011
                    position0.Pon((length(position0.Pon)+1):length(position.qrs))=NaN;
                    position0.P((length(position0.P)+1):length(position.qrs))=NaN;
                    position0.Poff((length(position0.Poff)+1):length(position.qrs))=NaN;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% finishing excerpt %%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Last annotated position
            if intervalo(2)>0
                if size(sig,2)>1
                    if ~isnan(position.Toff(min(intervalo(2),length(position.Toff)))),
                        ultimo_anot = position.Toff(min(intervalo(2),length(position.Toff)));
                        ultimo_anottyp=10;
                    elseif ~isnan(position.Tprima(min(intervalo(2),length(position.Tprima)))),
                        ultimo_anot = position.Tprima(min(intervalo(2),length(position.Tprima)));
                        ultimo_anottyp=9;
                    elseif ~isnan(position.T(min(intervalo(2),length(position.T)))),
                        ultimo_anot = position.T(min(intervalo(2),length(position.T)));
                        ultimo_anottyp=8;
                    elseif ~isnan(position.QRSoff(min(intervalo(2),length(position.QRSoff)))),
                        ultimo_anot = position.QRSoff(min(intervalo(2),length(position.QRSoff)));
                        ultimo_anottyp=6;
                    else
                        ultimo_anot = ceil(position.qrs(min(intervalo(2),length(position.qrs))));
                        ultimo_anottyp=5;
                    end
                else
                    if ~isnan(position1.Toff(min(intervalo(2),length(position1.Toff)))),
                        ultimo_anot = position1.Toff(min(intervalo(2),length(position1.Toff)));
                        ultimo_anottyp=10;
                    elseif ~isnan(position1.Tprima(min(intervalo(2),length(position1.Tprima)))),
                        ultimo_anot = position1.Tprima(min(intervalo(2),length(position1.Tprima)));
                        ultimo_anottyp=9;
                    elseif ~isnan(position1.T(min(intervalo(2),length(position1.T)))),
                        ultimo_anot = position1.T(min(intervalo(2),length(position1.T)));
                        ultimo_anottyp=8;
                    elseif ~isnan(position1.QRSoff(min(intervalo(2),length(position1.QRSoff)))),
                        ultimo_anot = position1.QRSoff(min(intervalo(2),length(position1.QRSoff)));
                        ultimo_anottyp=6;
                    else
                        ultimo_anot = ceil(position1.qrs(min(intervalo(2),length(position1.qrs))));
                        ultimo_anottyp=5;
                    end
                end
            end
        end %%%%%%%%%%%%%%%%introduzido a 17/05/02 Rute

    end
    
    
    inisamp = endsamp +2-endoverlap-begoverlap;
end

if( ~messages.setup.wavedet.QRS_detection_only )

    eval(['position' str '.Pon((length(position' str '.Pon)+1):length(position' str '.qrs))=NaN;']);
    eval(['position' str '.P((length(position' str '.P)+1):length(position' str '.qrs))=NaN;']);
    eval(['position' str '.Pprima((length(position' str '.Pprima)+1):length(position' str '.qrs))=NaN;']);
    eval(['position' str '.Poff((length(position' str '.Poff)+1):length(position' str '.qrs))=NaN;']);
    eval(['position' str '.Ton((length(position' str '.Ton)+1):length(position' str '.qrs))=NaN;']);
    eval(['position' str '.T((length(position' str '.T)+1):length(position' str '.qrs))=NaN;']);
    eval(['position' str '.Tprima((length(position' str '.Tprima)+1):length(position' str '.qrs))=NaN;']);
    eval(['position' str '.Toff((length(position' str '.Toff)+1):length(position' str '.qrs))=NaN;']);
    eval(['position' str '.Ttipo((length(position' str '.Ttipo)+1):length(position' str '.qrs))=NaN;']);
    eval(['position' str '.Tscale((length(position' str '.Tscale)+1):length(position' str '.qrs))=NaN;']);
    eval(['position' str '.Pon(position' str '.Pon>=(position' str '.qrs))=NaN;']);
    eval(['position' str '.P(position' str '.P>=(position' str '.qrs))=NaN;']);
    eval(['position' str '.Pprima(position' str '.Pprima>=(position' str '.qrs))=NaN;']);
    eval(['position' str '.Poff(position' str '.Poff>=(position' str '.qrs))=NaN;']);
    eval(['position' str '.Pon(position' str '.Pon>=(position' str '.QRSon))=NaN;']);
    eval(['position' str '.P(position' str '.P>=(position' str '.QRSon))=NaN;']);
    eval(['position' str '.Pprima(position' str '.Pprima>=(position' str '.QRSon))=NaN;']);
    eval(['position' str '.Poff(position' str '.Poff>=(position' str '.QRSon))=NaN;']);
    eval(['position' str '.Pon(position' str '.Pon>=(position' str '.P))=NaN;']);
    eval(['position' str '.Poff(position' str '.Poff<=(position' str '.P))=NaN;']);
    eval(['position' str '.Ton(position' str '.Ton>=(position' str '.T))=NaN;']);
    eval(['position' str '.Toff(position' str '.Toff<=(position' str '.T))=NaN;']);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% position assigment %%%%%%%%%%%%%
% Remove void annotations at the end
%strf = {'Pon';'P';'Pprima';'Ptipo';'Poff';'QRSon';'Q';'R';'qrs';'Rprima';'S';'QRSoff';...
%    'QRSpa';'QRSpp';'Ton';'T';'Tprima';'Toff';'Ttipo';'Tscale';'QRSmainpos';'QRSmaininv'};
strf = [];
for i = 1:3
    eval(['strf =fieldnames(position' num2str(i) ');']); % RUTE 4 AGO2011
    eval(['aux = find(position' num2str(i) '.qrs == 0| isnan(position' num2str(i) '.qrs));']);
    
    for j = 1:length(strf)
        eval(['aux_l=aux(aux<=size(position' num2str(i) '.' strf{j} ',2));'])
        eval(['position' num2str(i) '.' strf{j} '(aux_l) = [];']);
    end
end

if size(lead,2)==1
    position=position1; % 28FEB2012
    if ~isempty(position.QRSmainpos) && ~isempty(position.QRSmaininv)
        a = length(find(position.QRSmainpos == position.qrs));
        b = length(find(position.QRSmaininv == position.qrs));
        if a >= b
            position.qrs_hrv = '+';
        else
            position.qrs_hrv = '-';
        end
    end
else
    %     strf = {'Pon';'P';'Poff';'QRSon';'Q';'R';'qrs';'Rprima';'S';'QRSoff';...
    %         'Ton';'T';'Tprima';'Toff';'Ttipo';'Tscale'};
    %     strf2 = {'Ttipoon';'Ttipooff'};
    %     strf3 = {'QRSonsetcriteria';'QRSoffcriteria';'contadorToff';...
    %         'R_inQRSon';'R_inQRSoff'};
    strn = {'0';''};
    for i = 1:2
        strM = 'MM = max([';
        eval(['aux = find(position' strn{i} '.qrs == 0 | isnan(position' strn{i} '.qrs));']);
        eval(['strf =fieldnames(position' strn{i} ');']); % RUTE 4 AGO2011
        strd = ',';
        for j = 1:length(strf)
            eval(['aux_l=aux(aux<=size(position' strn{i} '.' strf{j} ',2));'])
            eval(['position' strn{i} '.' strf{j} '(aux_l) = [];']);
            if j == length(strf)
                strd = ']);';
            end
            strM = [strM 'size(position' strn{i} '.' strf{j} ')' strd]; %#ok<AGROW>
        end
        eval(strM);
        for j = 1:length(strf)
            eval(['position' strn{i} '.' strf{j} '(:,end+1:MM)=NaN*ones(size(position' strn{i} '.' strf{j} ',1),MM-size(position' strn{i} '.' strf{j} ',2));'])
        end
    end
end
% try
%     annstruct = pos2ann(position);
%     if ~isempty(annstruct.time)
%         writeannot([matdir ecgnr '.' anot], annstruct);   % mexfile, faster!
%     end
% catch me
%     me.message = 'No MIT format annotation file saved.';
%     messages.warnings=[messages.warnings {me.message}];
% end
for j = 1:4
    eval(['positionaux.position' num2str(j-1) ' = position' num2str(j-1) ';']); % RUTE 09.07.2007
end
% clear global OPT ecgnr
if exist('fid','var') && fid > 2
    fclose(fid);  % JBOLEA 11MAY2010
end