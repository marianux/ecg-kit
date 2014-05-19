% examples of wavedet_3D usage
%
%
%   Created by Rute Almeida (rbalmeid@unizar.es) on OCT, 2008.
%   Modified: NOV2011
%

% First define/correct the path for data
sigdir='E:\dados\QT\'; % path to the signal files
headir=sigdir; % path to the header files
matdir='E:\dados\QT\multilead_new\';  % path to the annotation files (to be saved)

%%%%% Examples
% QT database example (same as MIT-BIH)
registros =char('31','32');
ft=0; % format file: 0 for MIT header
anot=''; % inicial part of annotation files extension 
lead= 1; % lead - just one for single lead use
%lead=[1 2] % for multilead 2D lead
messages = new_empty_mesg;
messages.ss=[];
messages.setup.filter_bank_design=1;
for irec=1%:size(registros,1);
    ecgnr= [ 'sel' deblank(registros(irec,:))];
    %[position, aux, messages1] = wavedet_3D(sigdir,headir,matdir,ecgnr,ft,anot,lead);
    messages.setup.filter_bank_design=0;
    [position_n0, aux_n0, messages1_n0] = wavedet_3D(sigdir,headir,matdir,ecgnr,ft,anot,lead,[],[],[],[],[],messages);
    messages.setup.filter_bank_design=1;
    [position_n1, aux_n1, messages1_n1] = wavedet_3D(sigdir,headir,matdir,ecgnr,ft,anot,lead,[],[],[],[],[],messages);
end
marks1 = position_n1.qrs;
marks0 = position_n0.qrs;


% ESA example with external reading
setup.sigdir='E:\dados\ESA\AO-06-BR-ST-DLR\Mortara_Format\F1_BCD-5\'; %signal folder
setup.matdir='E:\dados\ESA\AO-06-BR-ST-DLR\Mortara_Format\F1_BCD-5\'; %output folder
setup.name='Hour1RawData'; %file to process
setup.ft=4; %data format
anot_label='';
inpt.heasig.freq=1000;
% for delinetion with internal QRS detection and no transformations
setup.flags=[0 0]; % [qrs_flag leadsynth_flag flagaux nsamp flag_sig_quality] 
lead{2}=1:12;% leads to use in the transformation
setup.t=[inpt.heasig.freq*9*60 inpt.heasig.freq*10*60]; % time to process in samples
setup.aname=[];
setup.dirann=[];
messages = new_empty_mesg;
messages.setup.wavedet.finvent_tol=0.15; % for high samplig frequencies eg stress test
if isunix, sep = '/'; else sep = '\'; end

% to  read a signal
[sig_all,messages,inpt.heasig] = readsignal(setup.sigdir,setup.name,[1 setup.t(2)],setup.ft,[],lead{2},messages);
inpt.heasig.d=1:inpt.heasig.nsig;
inpt.iBdir = [setup.matdir filesep];

for l=1:length(lead{2})
    inpt.ecgannot=num2str(lead{2}(l));
    inpt.annotmat = [setup.name '.' inpt.ecgannot]; % RUTE 06MAY2010
    lead{1}=lead{2}(l); %#ok<SAGROW> %lead to delineate:
    [position{l} aux messages1] = wavedet_3D('','', setup.matdir, setup.name,...
        setup.ft,anot_label, lead, setup.t, setup.flags, inpt.ecgannot,...
        setup.aname, setup.dirann,messages,sig_all,inpt.heasig); %#ok<SAGROW,NASGU>
    genera_results(inpt,position{l});
end

% SLR after SL
marks = pos2matrix(position);
[marksSLR ind messages1] = SLR(marks,inpt.heasig.freq);
