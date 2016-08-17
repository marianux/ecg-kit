%% Reads ECG recording in HL7a format
% Reads ECG recordings in Hl7a format from the Chinese database (CCDD). Implements the
% standard "HL7 aECG Implementation Guide - March 21, 2005" available in: 
% 
% https://www.hl7.org/documentcenter/public_temp_75706E59-1C23-BA17-0C25A0CC0545890C/wg/rcrim/annecg/aECG%20Implementation%20Guide%202005-03-21%20final%203.pdf
% 
% Arguments:
%   + filename: recording to be read.
%   + start_sample: (opt) start sample to read. Default 1.
%   + end_sample: (opt) end sample to read. Default min(All recording, ECG block of 200 Mbytes)
% 
% Output:
%   + ECG: the ECG block
%   + heasig: header with the ECG properties. 
%   + ann: annotations for the ECG recordings.
% 
% Limits:
% This routine is limited to read blocks smaller than 200 Mbytes for
% performance reasons. You can disable this limit by doing:
% MaxIOread = Inf; %megabytes
% 
% See also read_ishne_ann, read_ishne_header, read_ECG, ECGwrapper
% 
% Author: Mariano Llamedo Soria
% <matlab:web('mailto:llamedom@electron.frba.utn.edu.ar','-browser') (email)> 
% Version: 0.1 beta
% Birthdate: 02/05/2016
% Last update: 02/05/2016
% Copyright 2008-2016
% 
function [ ECG, heasig, ann, single_lead_positions, last_sample ] = read_hl7a_format(filename, start_sample, end_sample)

%% Tables

%% Lead Names
cLeadNamesHL7a = { ...
                    'MDC_ECG_LEAD_CONFIG', 'Unspecified lead', 'Unspecified lead'; ...
                    'MDC_ECG_LEAD_I', 'I', 'Lead I'; ...
                    'MDC_ECG_LEAD_II', 'II', 'Lead II'; ...
                    'MDC_ECG_LEAD_V1', 'V1', 'Lead V1'; ...
                    'MDC_ECG_LEAD_V2', 'V2', 'Lead V2'; ...
                    'MDC_ECG_LEAD_V3', 'V3', 'Lead V3'; ...
                    'MDC_ECG_LEAD_V4', 'V4', 'Lead V4'; ...
                    'MDC_ECG_LEAD_V5', 'V5', 'Lead V5'; ...
                    'MDC_ECG_LEAD_V6', 'V6', 'Lead V6'; ...
                    'MDC_ECG_LEAD_V7', 'V7', 'Lead V7'; ...
                    'MDC_ECG_LEAD_V2R', 'V2R', 'Lead V2Reversed'; ...
                    'MDC_ECG_LEAD_V3R', 'V3R', 'Lead V3Reversed'; ...
                    'MDC_ECG_LEAD_V4R', 'V4R', 'Lead V4Reversed'; ...
                    'MDC_ECG_LEAD_V5R', 'V5R', 'Lead V5Reversed'; ...
                    'MDC_ECG_LEAD_V6R', 'V6R', 'Lead V6Reversed'; ...
                    'MDC_ECG_LEAD_V7R', 'V7R', 'Lead V7Reversed'; ...
                    'MDC_ECG_LEAD_X', 'X', 'Frank´s X'; ...
                    'MDC_ECG_LEAD_Y', 'Y', 'Frank´s Y'; ...
                    'MDC_ECG_LEAD_Z', 'Z', 'Frank´s Z'; ...
                    'MDC_ECG_LEAD_CC5', 'CC5', 'CC5 per V5 and V5R placement'; ...
                    'MDC_ECG_LEAD_CM5', 'CM5', 'CM5 per V5 placement'; ...
                    'MDC_ECG_LEAD_LA', 'LA', 'Left Arm'; ...
                    'MDC_ECG_LEAD_RA', 'Lead RA', 'Right Arm'; ...
                    'MDC_ECG_LEAD_LL', 'LL', 'Left Leg'; ...
                    'MDC_ECG_LEAD_fI', 'I', 'I'; ...
                    'MDC_ECG_LEAD_fE', 'E', 'E'; ...
                    'MDC_ECG_LEAD_fC', 'C', 'C'; ...
                    'MDC_ECG_LEAD_fA', 'A', 'A'; ...
                    'MDC_ECG_LEAD_fM', 'M', 'M'; ...
                    'MDC_ECG_LEAD_fF', 'F', 'F'; ...
                    'MDC_ECG_LEAD_fH', 'H', 'H'; ...
                    'MDC_ECG_LEAD_III', 'III', 'III'; ...
                    'MDC_ECG_LEAD_AVR', 'aVR', 'aVR augmented voltage right'; ...
                    'MDC_ECG_LEAD_AVL', 'aVL', 'aVL augmented voltage left'; ...
                    'MDC_ECG_LEAD_AVF', 'aVF', 'aVF augmented voltage foot'; ...
                    'MDC_ECG_LEAD_AVRneg', '–aVR', '?aVR'; ...
                    'MDC_ECG_LEAD_V8', 'V8', 'V8'; ...
                    'MDC_ECG_LEAD_V9', 'V9', 'V9'; ...
                    'MDC_ECG_LEAD_V8R', 'V8R', 'V8R'; ...
                    'MDC_ECG_LEAD_V9R', 'V9R', 'V9R'; ...
                    'MDC_ECG_LEAD_D', 'D', 'D (Nehb – Dorsal)'; ...
                    'MDC_ECG_LEAD_A', 'A', 'A (Nehb – Anterior)'; ...
                    'MDC_ECG_LEAD_J', 'J', 'J (Nehb – Inferior)'; ...
                    'MDC_ECG_LEAD_DEFIB', 'Defib', 'Defibrillator lead: anterior-lateral'; ...
                    'MDC_ECG_LEAD_EXTERN', 'Extern', 'External pacing lead: anterior-posterior'; ...
                    'MDC_ECG_LEAD_A1', 'A1', 'A1 (Auxiliary unipolar lead #1)'; ...
                    'MDC_ECG_LEAD_A2', 'A2', 'A2 (Auxiliary unipolar lead #2)'; ...
                    'MDC_ECG_LEAD_A3', 'A3', 'A3 (Auxiliary unipolar lead #3)'; ...
                    'MDC_ECG_LEAD_A4', 'A4', 'A4 (Auxiliary unipolar lead #4)'; ...
                    'MDC_ECG_LEAD_C', 'Chest', 'Chest lead'; ...
                    'MDC_ECG_LEAD_V', 'V', 'Precordial lead'; ...
                    'MDC_ECG_LEAD_VR', 'VR', 'VR nonaugmented voltage vector of RA'; ...
                    'MDC_ECG_LEAD_VL', 'VL', 'VL nonaugmented voltage vector of LA'; ...
                    'MDC_ECG_LEAD_VF', 'VF', 'VF nonaugmented voltage vector of LL'; ...
                    'MDC_ECG_LEAD_MCL', 'MCL', 'Modified chest lead (left arm indifferent)'; ...
                    'MDC_ECG_LEAD_MCL1', 'MCL1', 'MCL per V1 placement'; ...
                    'MDC_ECG_LEAD_MCL2', 'MCL2', 'MCL per V2 placement'; ...
                    'MDC_ECG_LEAD_MCL3', 'MCL3', 'MCL per V3 placement'; ...
                    'MDC_ECG_LEAD_MCL4', 'MCL4', 'MCL per V4 placement'; ...
                    'MDC_ECG_LEAD_MCL5', 'MCL5', 'MCL per V5 placement'; ...
                    'MDC_ECG_LEAD_MCL6', 'MCL6', 'MCL per V6 placement'; ...
                    'MDC_ECG_LEAD_CC', 'CC', 'Chest lead (symmetric placement)'; ...
                    'MDC_ECG_LEAD_CC1', 'CC1', 'CC1 per V1 and V1R placement'; ...
                    'MDC_ECG_LEAD_CC2', 'CC2', 'CC2 per V2 and V2R placement'; ...
                    'MDC_ECG_LEAD_CC3', 'CC3', 'CC3 per V3 and V3R placement'; ...
                    'MDC_ECG_LEAD_CC4', 'CC4', 'CC4 per V4 and V4R placement'; ...
                    'MDC_ECG_LEAD_CC6', 'CC6', 'CC6 per V6 and V6R placement'; ...
                    'MDC_ECG_LEAD_CC7', 'CC7', 'CC7 per V7 and V8R placement'; ...
                    'MDC_ECG_LEAD_CM', 'CM', 'Chest-manubrium'; ...
                    'MDC_ECG_LEAD_CM1', 'CM1', 'CM1 per V1 placement'; ...
                    'MDC_ECG_LEAD_CM2', 'CM2', 'CM2 per V2 placement'; ...
                    'MDC_ECG_LEAD_CM3', 'CM3', 'CM3 per V3 placement'; ...
                    'MDC_ECG_LEAD_CM4', 'CM4', 'CM4 per V4 placement'; ...
                    'MDC_ECG_LEAD_CM6', 'CM6', 'CM6 per V6 placement'; ...
                    'MDC_ECG_LEAD_CM7', 'CM7', 'CM7 per V7 placement'; ...
                    'MDC_ECG_LEAD_CH5', 'CH5', '-'; ...
                    'MDC_ECG_LEAD_CS5', 'CS5', 'negative: right infraclavicular fossa'; ...
                    'MDC_ECG_LEAD_CB5', 'CB5', 'negative: low right scapula'; ...
                    'MDC_ECG_LEAD_CR5', 'CR5', '-'; ...
                    'MDC_ECG_LEAD_ML', 'ML', 'ML modified limb lead ~ Lead II'; ...
                    'MDC_ECG_LEAD_AB1', 'AB1', 'AB1 (auxiliary bipolar lead #1)'; ...
                    'MDC_ECG_LEAD_AB2', 'AB2', 'AB2 (auxiliary bipolar lead #2)'; ...
                    'MDC_ECG_LEAD_AB3', 'AB3', 'AB3 (auxiliary bipolar lead #3)'; ...
                    'MDC_ECG_LEAD_AB4', 'AB4', 'AB4 (auxiliary bipolar lead #4)'; ...
                    'MDC_ECG_LEAD_ES', 'ES', 'EASI™ ES'; ...
                    'MDC_ECG_LEAD_AS', 'AS', 'EASI AS'; ...
                    'MDC_ECG_LEAD_AI', 'AI', 'EASI AI'; ...
                    'MDC_ECG_LEAD_S', 'S', 'EASI upper sternum lead'; ...
                    'MDC_ECG_LEAD_dI', 'dI', 'derived lead I'; ...
                    'MDC_ECG_LEAD_dII', 'dII', 'derived lead II'; ...
                    'MDC_ECG_LEAD_dIII', 'dIII', 'derived lead III'; ...
                    'MDC_ECG_LEAD_daVR', 'daVR', 'derived lead aVR'; ...
                    'MDC_ECG_LEAD_daVL', 'daVL', 'derived lead aVL'; ...
                    'MDC_ECG_LEAD_daVF', 'daVF', 'derived lead aVF'; ...
                    'MDC_ECG_LEAD_dV1', 'dV1', 'derived lead V1'; ...
                    'MDC_ECG_LEAD_dV2', 'dV2', 'derived lead V2'; ...
                    'MDC_ECG_LEAD_dV3', 'dV3', 'derived lead V3'; ...
                    'MDC_ECG_LEAD_dV4', 'dV4', 'derived lead V4'; ...
                    'MDC_ECG_LEAD_dV5', 'dV5', 'derived lead V5'; ...
                    'MDC_ECG_LEAD_dV6', 'dV6', 'derived lead V6'; ...
                    'MDC_ECG_LEAD_RL', 'RL', 'right leg'; ...
                    'MDC_ECG_LEAD_CV5RL', 'CV5RL', 'Canine fifth right intercostals space near the edge of the sternum at the most curved part of the costal cartilage'; ...
                    'MDC_ECG_LEAD_CV6LL', 'CV6LL', 'Canine sixth left intercostals space near the edge of the sternum at the most curved part of the costal cartilage'; ...
                    'MDC_ECG_LEAD_CV6LU', 'CV6LU', 'Canine sixth left intercostals space at the costochondral junction'; ...
                    'MDC_ECG_LEAD_V10', 'V10', 'Canine over dorsal spinous process of the seventh thoracic vertebra'; ...
                    };

%% Wave Names

cWaveNamesHL7a = { ...
                    'MDC_ECG_WAVC_PWAVE', 'P', 'P wave'; ...
                    'MDC_ECG_WAVC_PPWAVE', 'P´', 'P´ wave (second deflection in P wave) (P and P´ waves have opposite signs)'; ...
                    'MDC_ECG_WAVC_PPPWAVE', 'P"', 'P" wave (third deflection in P wave) (P´ and P" waves have opposite signs)'; ...
                    'MDC_ECG_WAVC_QWAVE', 'Q', 'Q wave'; ...
                    'MDC_ECG_WAVC_QSWAVE', 'QS', 'QS wave'; ...
                    'MDC_ECG_WAVC_RWAVE', 'R', 'R wave'; ...
                    'MDC_ECG_WAVC_RRWAVE', 'R´', 'R´ wave (second deflection in R Wave) (R and R´ have same sign)'; ...
                    'MDC_ECG_WAVC_RRRWAVE', 'R"', 'R" wave (third deflection in R Wave) (R, R´ and R" have same sign)'; ...
                    'MDC_ECG_WAVC_NOTCH', 'Notch', 'Notch a slight but distinct change in the direction of a WAVC deflection, contained entirely within that deflection. Typically associated with Q-, R- and/or S-wave.'; ...
                    'MDC_ECG_WAVC_SWAVE', 'S', 'S wave (S and R/R´ waves have opposite signs)'; ...
                    'MDC_ECG_WAVC_SSWAVE', 'S´', 'S´ wave (second deflection in S Wave) (S´ and R´/R" waves have opposite signs)'; ...
                    'MDC_ECG_WAVC_SSSWAVE', 'S"', 'S" wave (third deflection in S Wave) (S´ and R´/R" waves have opposite signs)'; ...
                    'MDC_ECG_WAVC_TWAVE', 'T', 'T wave'; ...
                    'MDC_ECG_WAVC_TTWAVE', 'T´', 'T´ wave (second deflection in T Wave) (T and T´ waves have opposite signs)'; ...
                    'MDC_ECG_WAVC_UWAVE', 'U', 'U wave'; ...
                    'MDC_ECG_WAVC_DELTA', 'Delta', 'Delta wave'; ...
                    'MDC_ECG_WAVC_IWAVE', 'I', 'Isoelectric region between global QRS onset and actual onset of QRS in given lead'; ...
                    'MDC_ECG_WAVC_KWAVE', 'K', 'Isoelectric region between actual offset of QRS in given lead and global QRS offset'; ...
                    'MDC_ECG_WAVC_JWAVE', 'J', 'Osborne wave, late and typically upright terminal deflection of QRS complex; amplitude increases as temperature declines. ECG finding typically associated with hypothermia.'; ...
                    'MDC_ECG_WAVC_PQRSTWAVE', 'PQRST', 'Entire Beat (Pon to Toff, excluding U)'; ...
                    'MDC_ECG_WAVC_QRSTWAVE', 'QRST', 'Entire Beat (Qon to Toff, excluding P and U)'; ...
                    'MDC_ECG_WAVC_QRSWAVE', 'QRS', 'Entire QRS (excluding P, T and U)'; ...
                    'MDC_ECG_WAVC_TUWAVE', 'TU', 'TU fused wave'; ...
                    'MDC_ECG_WAVC_VFLWAVE', 'V flutter', 'wave Ventricular flutter wave (optional) (the appropriate ventricular rhythm call is mandatory)'; ...
                    'MDC_ECG_WAVC_AFLWAVE', 'Atrial flutter', 'wave Atrial flutter wave (optional) (the appropriate atrial rhythm call is mandatory)'; ...
                    'MDC_ECG_WAVC_ISO', 'Isoelectric point', 'Isoelectric point or segment'; ...
                    'MDC_ECG_WAVC_PRSEG', 'PR Segment', 'PR Segment'; ...
                    'MDC_ECG_WAVC_STSEG', 'ST Segment', 'ST Segment'; ...
                    'MDC_ECG_WAVC_STJ', 'J-point', 'J-point'; ...
                    'MDC_ECG_WAVC_STM', 'ST meas', 'point ST measurement point'; ...
                    'MDC_ECG_WAVC_ARFCT', 'Artifact Isolated', 'qrs-like artifact'; ...
                    'MDC_ECG_WAVC_CALP', 'Calibration pulse', 'Calibration pulse (individual pulse)'; ...
                    'MDC_ECG_WAVC_STCH', 'ST change', 'ST change'; ...
                    'MDC_ECG_WAVC_TCH', 'T-wave change', 'T-wave change'; ...
                    'MDC_ECG_WAVC_VAT', 'Ventricular Activation', 'Time Ventricular Activation Time also termed the intrinsic (or intrinsicoid) deflection onset to peak of depolarization wave.'; ...
                    };

ecgkit_wave_defs;

[~, wave_subset_idx, aux_idx] = intersect(cWaveNamesHL7a(:,1), cHL7aECG_translation(:,1) );
cHL7aECG_translation = cHL7aECG_translation(aux_idx,:);
[~, QRwaves_idx] = intersect(cWaveNamesHL7a(wave_subset_idx,1), {'MDC_ECG_WAVC_QWAVE', 'MDC_ECG_WAVC_RWAVE', 'MDC_ECG_WAVC_QSWAVE', 'MDC_ECG_WAVC_QRSWAVE'} );
                
%% Beat types

% The second column is an arbitrary asignation of HL7a beats to AAMI EC-57
% standard.
cBeatTypesHL7a = { ...
                    'MDC_ECG_BEAT', 'Q', 'Any beat (unspecified; included in heart rate)'; ...
                    'MDC_ECG_BEAT_NORMAL', 'N', 'Normal Beat Normal beat (sinus beat normal conduction)'; ...
                    'MDC_ECG_BEAT_ABNORMAL', 'Q', 'Abnormal Beat Abnormal beat'; ...
                    'MDC_ECG_BEAT_DOMINANT', 'Q', 'Dominant Beat Dominant beat (typically normal but may not be) (predominant morphology typically used for ST measurement)'; ...
                    'MDC_ECG_BEAT_SV_P_C', 'S', 'Supraventricular premature contraction Supraventricular premature contraction (atrial or nodal premature beat with normal QRS morphology)'; ...
                    'MDC_ECG_BEAT_ATR_P_C', 'S', 'Atrial Premature contraction Atrial premature contraction (beat)'; ...
                    'MDC_ECG_BEAT_JUNC_P_C', 'S', 'Junctional premature contraction Junctional (nodal) premature contraction'; ...
                    'MDC_ECG_BEAT_ATR_P_C_ABER', 'S', 'Aberrated atrial premature beat Aberrated atrial premature beat (Ashman beat) (atrial premature beat with abnormal QRS morphology)'; ...
                    'MDC_ECG_BEAT_R', 'S', 'Aberrated atrial premature beat Aberrated atrial premature beat (Ashman beat) (atrial premature beat with abnormal QRS morphology)'; ...
                    'MDC_ECG_BEAT_ATR_PWAVE_B', 'Q', 'Non-conducted p-wave Non-conducted p-wave (blocked)'; ...
                    'MDC_ECG_BEAT_LK', 'Q', 'Non-conducted p-wave Non-conducted p-wave (blocked)'; ...
                    'MDC_ECG_BEAT_V_P_C', 'V', 'Ventricular premature contraction Ventricular premature contraction (beat)'; ...
                    'MDC_ECG_BEAT_V_P_C_FUSION', 'F', 'Fusion of ventricular and normal beat Fusion of ventricular and normal beat'; ...
                    'MDC_ECG_BEAT_V_P_C_RonT', 'V', 'R-on-T premature ventricular beat R-on-T premature ventricular beat'; ...
                    'MDC_ECG_BEAT_SV_ESC', 'S', 'Supraventricular escape beat (least specific)'; ...
                    'MDC_ECG_BEAT_ATR_ESC', 'S', 'Atrial escape beat'; ...
                    'MDC_ECG_BEAT_JUNC_ESC', 'S', 'Junctional (nodal) escape beat'; ...
                    'MDC_ECG_BEAT_V_ESC', 'V', 'Ventricular escape beat'; ...
                    'MDC_ECG_BEAT_BB_BLK', 'N', 'bundle branch block beat (unspecified)'; ...
                    'MDC_ECG_BEAT_LBB_BLK_COMP', 'N', 'left bundle branch block beat'; ...
                    'MDC_ECG_BEAT_LBB_BLK_INCOMP', 'N', 'incomplete left bundle branch block beat'; ...
                    'MDC_ECG_BEAT_RBB_BLK_COMP', 'N', 'right bundle branch block beat'; ...
                    'MDC_ECG_BEAT_RBB_BLK_INCOMP', 'N', 'incomplete right bundle branch block beat'; ...
                    'MDC_ECG_BEAT_BLK_ANT_L_HEMI', 'N', 'left anterior fascicular block beat (common)'; ...
                    'MDC_ECG_BEAT_BLK_POS_L_HEMI', 'N', 'left posterior fascicular block beat (rare)'; ...
                    'MDC_ECG_BEAT_BLK_BIFASC', 'N', 'bifascicular block beat'; ...
                    'MDC_ECG_BEAT_BLK_TRIFASC', 'N', 'trifascicular block beat'; ...
                    'MDC_ECG_BEAT_BLK_BILAT', 'N', 'bilateral bundle-branch block beat'; ...
                    'MDC_ECG_BEAT_BLK_IVCD', 'N', 'intraventricular conduction disturbance (non-specific block)'; ...
                    'MDC_ECG_BEAT_PREX', 'Q', 'pre-excitation (least specific)'; ...
                    'MDC_ECG_BEAT_WPW_UNK', 'Q', 'Wolf-Parkinson-White syndrome (less specific)'; ...
                    'MDC_ECG_BEAT_WPW_A', 'Q', 'Wolf-Parkinson type A'; ...
                    'MDC_ECG_BEAT_WPW_B', 'Q', 'Wolf-Parkinson type B'; ...
                    'MDC_ECG_BEAT_LGL', 'Q', 'Lown-Ganong-Levine syndrome'; ...
                    'MDC_ECG_BEAT_PACED', 'Q', 'Paced beat (with ventricular capture)'; ...
                    'MDC_ECG_BEAT_PACED_FUS', 'Q', 'Pacemaker Fusion beat'; ...
                    'MDC_ECG_BEAT_UNKNOWN', 'Q', 'Unclassifiable beat'; ...
                    'MDC_ECG_BEAT_LEARN', 'Q', 'Learning (beat during initial learning phase)'; ...
                    };

                
%% Code
%No leer bloques mas grandes de 200 megabytes
MaxIOread = 200; %megabytes

if( nargin < 2 || isempty( start_sample ) )
    start_sample = 1;
else
    start_sample = max(1,start_sample);
end

ann = [];
heasig = [];
ECG = [];
last_sample = [];

if( nargout > 1 )
    bHeaderRequired = true;
else
    bHeaderRequired = false;
end

if( nargout > 2 )
    bAnnRequired = true;
else
    bAnnRequired = false;
end

try
   xDoc = xmlread(filename);
catch
    error('read_hl7a_format:ReadError', 'Failed to read XML file %s.\n', filename);
end

% get series in the file
allSeries = xDoc.getElementsByTagName('series');

if(allSeries.getLength == 0)
    error('read_hl7a_format:NoSeries', 'No series found in %s.\n', filename);
elseif(allSeries.getLength > 1)
    disp_string_framed(2, sprintf( 'More than one serie in %s. Reading only the first one', filename) );
end

allSeries = allSeries.item(0);

if(bHeaderRequired)

    [~, heasig.recname ] = fileparts(filename);

    efectivetime = allSeries.getElementsByTagName('effectiveTime');
    efectivetime = efectivetime.item(0);

    loww = efectivetime.getElementsByTagName('low');
    low_val = get_time_from_xml_tag(loww.item(0), false);
    
%     loww = efectivetime.getElementsByTagName('high');
%     high_val = get_time_from_xml_tag(loww.item(0), false);
    
    if( isempty(low_val) )
        base_time = datenum('20000101000000.000', 'yyyymmddHHMMSS.FFF');
        % unknown generic date
        heasig.btime = '00:00:00';
        heasig.bdate = '01/01/2000';
    else
        base_time = low_val;
        heasig.bdate = datestr(low_val, 'yyyy/mm/dd');
        heasig.btime = datestr(low_val, 'HH:MM:SS');
    end
    
end
    
%% Signal Parsing

sequenceSet = allSeries.getElementsByTagName('sequenceSet');

node_index_limit = realmax;

if( sequenceSet.getLength > 1 )
    disp_string_framed(2, sprintf( 'More than one signal in %s. Reading only the first one', filename) );
    % several signals in this file, get the limit to stop searching data.
    aux_val = sequenceSet.item(1);
    node_index_limit = aux_val.getNodeIndex;
end

sequenceSet = sequenceSet.item(0);

allComponents = sequenceSet.getElementsByTagName('component');

heasig.nsig = allComponents.getLength-1;
heasig.desc = repmat({''}, heasig.nsig,1);
heasig.adczero = zeros(heasig.nsig,1);
heasig.gain = ones(heasig.nsig,1);
heasig.units = repmat({''}, heasig.nsig,1);

lead_idx = 1;
lead_translate_idx = zeros(size(cLeadNamesHL7a,1),1);

for ii = 0:(allComponents.getLength-1)

    thisComp = allComponents.item(ii);
    
    allVals = thisComp.getElementsByTagName('value');
    
    thisVal = allVals.item(0);
    
    thisVal_att = thisVal.getAttributes;
    
    for jj = 0:(thisVal_att.getLength-1)
       
        this_att = thisVal_att.item(jj);

        aux_val = char(this_att.getValue);
        
        if( strcmpi(aux_val, 'SLIST_PQ') )
        %% ECG leads
            
            if(bHeaderRequired)
        
                %% ADC level
                thisOrig = thisVal.getElementsByTagName('origin');
                thisOrig = thisOrig.item(0);

                thisOrig_att = thisOrig.getAttributes;

                for kk = 0:(thisOrig_att.getLength-1)

                    this_att = thisOrig_att.item(kk);
                    heasig.adczero(lead_idx) = 1;

                    aux_val = this_att.getName;

                    if( strcmpi(aux_val, 'value') )
                        heasig.adczero(lead_idx) = heasig.adczero(lead_idx) * str2double(char(this_att.getValue));

                    elseif( strcmpi(aux_val, 'unit') )

                        aux_val = char(this_att.getValue);

                        switch(aux_val)
                            case 'V'
                                heasig.adczero(lead_idx) = heasig.adczero(lead_idx) ;
                            case 'mV'
                                heasig.adczero(lead_idx) = heasig.adczero(lead_idx) * 1e-3;
                            case 'uV'
                                heasig.adczero(lead_idx) = heasig.adczero(lead_idx) * 1e-6;
                            otherwise
                                error('read_hl7a_format:ParseError', 'Parse error at %s. Check lead %d <origin unit = %s\n', filename, lead_idx, aux_val);
                        end
                    else
                        error('read_hl7a_format:ParseError', 'Parse error at %s. Check lead %d, attribute %s:%s\n', filename, lead_idx, char(this_att.getName), char(this_att.getValue));
                    end

                end            

                %% ADC Gain

                thisScale = thisVal.getElementsByTagName('scale');
                thisScale = thisScale.item(0);

                thisScale_att = thisScale.getAttributes;        

                thisScale_att = thisScale.getAttributes;

                for kk = 0:(thisScale_att.getLength-1)

                    this_att = thisScale_att.item(kk);
                    heasig.gain(lead_idx) = 1;

                    aux_val = this_att.getName;

                    if( strcmpi(aux_val, 'value') )
                        heasig.gain(lead_idx) = heasig.gain(lead_idx) * 1/str2double(this_att.getValue);

                    elseif( strcmpi(aux_val, 'unit') )

                        aux_val = char(this_att.getValue);
                        heasig.units{lead_idx} = aux_val;

    %                     switch(aux_val)
    %                         case 'V'
    %                             heasig.gain(lead_idx) = heasig.gain(lead_idx) ;
    %                         case 'mV'
    %                             heasig.gain(lead_idx) = heasig.gain(lead_idx) * 1e-3;
    %                         case 'uV'
    %                             heasig.gain(lead_idx) = heasig.gain(lead_idx) * 1e-6;
    %                         case 'nV'
    %                             heasig.gain(lead_idx) = heasig.gain(lead_idx) * 1e-9;
    %                         otherwise
    %                             error('read_hl7a_format:ParseError', 'Parse error at %s. Check lead %d <scale unit = %s\n', filename, lead_idx, aux_val);
    %                     end
                    else
                        error('read_hl7a_format:ParseError', 'Parse error at %s. Check lead %d, attribute %s:%s\n', filename, lead_idx, char(this_att.getName), char(this_att.getValue));
                    end

                end                     

    %             <origin value="0" unit="uV" />
    %             <scale value="4.76837" unit="uV" />
    %             digits

                %% Lead name

                thisCode = thisComp.getElementsByTagName('code');
                thisCode = thisCode.item(0);

                thisCode_att = thisCode.getAttributes;        

                for kk = 0:(thisCode_att.getLength-1)

                    this_att = thisCode_att.item(kk);

                    aux_val = this_att.getName;

                    if( strcmpi(aux_val, 'code') )
                        [aux_val, aux_idx ]= intersect(upper(cLeadNamesHL7a(:,1)), upper(char(this_att.getValue)));

                        if( isempty(aux_idx) )
                            error('read_hl7a_format:ParseError', 'Unknown lead name at %s. Check lead %s.\n', filename, char(this_att.getValue));
                        else
                            heasig.desc(lead_idx) = cLeadNamesHL7a(aux_idx,2);
                            lead_translate_idx(aux_idx) = lead_idx;
                        end
                    end

                end                     

            end
            %% Data samples
            
            thisDigits = thisVal.getElementsByTagName('digits');
            thisDigits = thisDigits.item(0);            

            if( isempty(ECG) )
                aux_val = char(thisDigits.getTextContent);
                aux_val = regexprep(aux_val, '[^-+\d]+', ' ');
                ECG = colvec(str2num(aux_val));
                heasig.nsamp = size(ECG,1);
                
                if( nargin < 3 || isempty( end_sample ) )
                    %Intento la lectura total por defecto
                    samples2read = heasig.nsamp - (start_sample-1);
                else
                    samples2read = min(heasig.nsamp, end_sample) - (start_sample-1);
                end
                
                end_sample = start_sample + samples2read - 1;
                heasig.nsamp = samples2read;              
                ECG = colvec(ECG(start_sample:end_sample));
                
            else
                aux_val = char(thisDigits.getTextContent);
                aux_val = colvec(str2num(regexprep(aux_val, '[^-+\d]+', ' ')));
                ECG(:,lead_idx) = colvec(aux_val(start_sample:end_sample));
            end
            
            lead_idx = lead_idx + 1;
            
        elseif( strcmpi(aux_val, 'GLIST_TS') )
        %% time sequence -> sampling rate
            if(bHeaderRequired)
        
                thisInc = thisVal.getElementsByTagName('increment');
                thisInc = thisInc.item(0);

                thisInc_att = thisInc.getAttributes;

                for kk = 0:(thisInc_att.getLength-1)

                    this_att = thisInc_att.item(kk);
                    heasig.freq = 1;

                    aux_val = char(this_att.getName);

                    if( strcmpi(aux_val, 'value') )
                        heasig.freq = 1/str2double(this_att.getValue);
                    elseif( strcmpi(aux_val, 'unit') )
                        aux_val = char(this_att.getValue);
                        switch(aux_val)
                            case 's'
                                heasig.freq = heasig.freq;
                            case 'ms'
                                heasig.freq = heasig.freq * 1e3;
                            case 'us'
                                heasig.freq = heasig.freq * 1e6;
                            otherwise
                                error('read_hl7a_format:ParseError', 'Parse error at %s. Check unit = %s\n', filename, aux_val);
                        end
                    else
                        error('read_hl7a_format:ParseError', 'Parse error at %s. Check %s:%s\n', filename, aux_val, char(this_att.getValue));
                    end

                end
            end
        end
        
    end
    
end

heasig.desc = char(heasig.desc);
heasig.units = char(heasig.units);

%% Annotations parsing

if(bAnnRequired)

    ann.time = [];
    ann.anntyp = [];
    
    bDesktop = usejava('desktop');
    
    if(bDesktop)            
        % Activate the progress_struct bar.
        pb = progress_bar(heasig.recname);
    end
    
    % empty output struct
    for fn = cAnnotationOutputFields
        single_lead_positions.(fn{1}) = [];
    end
    single_lead_positions = repmat(single_lead_positions, heasig.nsig, 1);
    
    % get annotationSet tag in the file
    annSets = xDoc.getElementsByTagName('annotationSet');
    
    next_index_base = 0;
    annSet_idx = 0;
    thisAnnSet = annSets.item(annSet_idx);
    
    bAnnSetFound = false;
    
    while( annSet_idx < annSets.getLength )

        thisAnnSet = annSets.item(annSet_idx);
        
        if( thisAnnSet.getNodeIndex > node_index_limit ) 
            break;
        end
        
        strAnnSet = sprintf( 'Annotation set %d', annSet_idx );
            
        AnnSetHash = thisAnnSet.hashCode;
        
        comp_idx = 0;
        allComponents = thisAnnSet.getElementsByTagName('component');

        pb.Loops2Do = max(1, allComponents.getLength-1);
        
        while( comp_idx < allComponents.getLength )

            if(bDesktop)
                %start of the progress_struct loop 0%
                pb.start_loop();
            end
            
            bCompParsedOk = false;
            
            thisComp = allComponents.item(comp_idx);

            if( thisComp.getNodeIndex > next_index_base )

                thisAnn = thisComp.getElementsByTagName('annotation');

                if( thisAnn.getLength > 0 )

                    thisAnn = thisAnn.item(0);

                    thisCode = thisAnn.getElementsByTagName('code');

                    thisValue = thisAnn.getElementsByTagName('value');

                    pb.checkpoint();

                    if( ( xml_tag_value(thisCode.item(0), 'code', 'MDC_ECG_WAVC' ) || xml_tag_value(thisCode.item(0), 'code', 'MDC_ECG_WAVC_TYPE' ) ) && ...
                            xml_tag_value(thisValue.item(1), 'code', 'MDC_ECG_WAVC_PEAK') ) 
                        % peak annotation

                        [~, strAux] = xml_tag_value(thisValue.item(0), 'code');

                        [~, wave_idx ] = intersect(cWaveNamesHL7a(wave_subset_idx,1), strAux);

                        if( ~isempty(wave_idx) )

                            % Annotation is global, unless specified by a "code"
                            % tag
                            bTimeRelative = true;
                            bAllleads = true;
                            code_idx = 2;
                            while(code_idx < thisCode.getLength)
                                [~, strAux] = xml_tag_value(thisCode.item(code_idx), 'code');

                                if( strcmpi(strAux, 'TIME_ABSOLUTE') )

                                    bTimeRelative = false;

                                else
                                    [~, lead_idx ]= intersect(upper(cLeadNamesHL7a(:,1)), upper(strAux));

                                    if( ~isempty(lead_idx) )
                                        bAllleads = false;
                                        break;
                                    end
                                end

                                code_idx = code_idx + 1;
                            end

                            this_peak = get_time_from_xml_tag(thisValue.item(2), bTimeRelative);

                            if( ~bTimeRelative )
                                % make annotations relative to the
                                % origin
                                this_peak = etime(datevec(this_peak), datevec(base_time));
                            end

                            if( bAllleads )
                                aux_idx = 1:heasig.nsig;
                            else
                                aux_idx = lead_translate_idx(lead_idx);
                            end

                            aux_val = round( this_peak * heasig.freq);

                            for ii = aux_idx

                                bCompParsedOk = true;

                                % peak
                                aux_idx = single_lead_positions(ii).(cHL7aECG_translation{wave_idx,2});
                                aux_idx = [aux_idx; aux_val];
                                single_lead_positions(ii).(cHL7aECG_translation{wave_idx,2}) = aux_idx;
                                bAnnSetFound = true;
                            end

                            aux_val = thisValue.item(2);
                            next_index_base = aux_val.getNodeIndex;

                        end


                    elseif( (   xml_tag_value(thisCode.item(0), 'code', 'MDC_ECG_BEAT' ) && ...
                            xml_tag_value(thisCode.item(1), 'code', 'MDC_ECG_WAVC' ) ...
                        ) || ...
                            xml_tag_value(thisCode.item(0), 'code', 'MDC_ECG_WAVC' )) 
                        % wave annotation (onset-offset)

                        [~, strAux ] = xml_tag_value(thisValue.item(0), 'code' );
                        [~, wave_idx ] = intersect(cWaveNamesHL7a(wave_subset_idx,1), strAux);
                        code_idx = 1;

                        if( isempty(wave_idx) )
                            code_idx = 2;
                            [~, strAux ] = xml_tag_value(thisValue.item(1), 'code' );
                            [~, wave_idx ] = intersect(cWaveNamesHL7a(wave_subset_idx,1), strAux);
                        end

                        if( ~isempty(wave_idx) )

                            % Annotation is global, unless specified by a "code"
                            % tag
                            bAllleads = true;
                            bTimeRelative = true;
                            while(code_idx < thisCode.getLength)
                                [~, strAux] = xml_tag_value(thisCode.item(code_idx), 'code');

                                if( strcmpi(strAux, 'TIME_ABSOLUTE') )

                                    bTimeRelative = false;

                                else

                                    [~, lead_idx ]= intersect(upper(cLeadNamesHL7a(:,1)), upper(strAux));

                                    if( ~isempty(lead_idx) )
                                        bAllleads = false;
                                        break;
                                    end

                                end

                                code_idx = code_idx + 1;
                            end

                            thisOnset = thisAnn.getElementsByTagName('low');
                            thisOffset = thisAnn.getElementsByTagName('high');

                            this_onset = get_time_from_xml_tag(thisOnset.item(0), bTimeRelative);
                            this_offset = get_time_from_xml_tag(thisOffset.item(0), bTimeRelative);

                            if( ~bTimeRelative )
                                % make annotations relative to the
                                % origin
                                if( ~isempty(this_onset) )
                                    this_onset = etime(datevec(this_onset), datevec(base_time));
                                end
                                if( ~isempty(this_offset) )
                                    this_offset = etime(datevec(this_offset), datevec(base_time));
                                end
                            end


                            if( bAllleads )
                                lead_idx2 = 1:heasig.nsig;
                            else
                                lead_idx2 = lead_translate_idx(lead_idx);
                            end

                            for ii = lead_idx2

                                bCompParsedOk = true;

                                if( ~isempty(cHL7aECG_translation{wave_idx,3}) )
                                    % onset
                                    aux_idx = single_lead_positions(ii).(cHL7aECG_translation{wave_idx,3});
                                    aux_val = round( this_onset * heasig.freq);
                                    aux_idx = [aux_idx; aux_val];
                                    single_lead_positions(ii).(cHL7aECG_translation{wave_idx,3}) = aux_idx;
                                    bAnnSetFound = true;
                                end
                                if( ~isempty(cHL7aECG_translation{wave_idx,4}) )
                                    % offset
                                    aux_idx = single_lead_positions(ii).(cHL7aECG_translation{wave_idx,4});
                                    aux_val = [aux_val; round( this_offset * heasig.freq)];
                                    wave_midpoint = round(mean(aux_val));
                                    aux_idx = [aux_idx; aux_val(end)];
                                    single_lead_positions(ii).(cHL7aECG_translation{wave_idx,4}) = aux_idx;
                                    bAnnSetFound = true;
                                end
                            end

                            if( any(QRwaves_idx == wave_idx) )
                                % heartbeat type
                                ann.time = [ann.time; wave_midpoint];
                                [~, strAux ] = xml_tag_value(thisValue.item(0), 'code' );
                                [~, beat_type_idx ]= intersect(cBeatTypesHL7a(:,1), upper(strAux));
                                if( isempty(beat_type_idx) )
                                    % Unclassified beat
                                    ann.anntyp = [ann.anntyp; cBeatTypesHL7a{1,2}];
                                else
                                    ann.anntyp = [ann.anntyp; cBeatTypesHL7a{beat_type_idx,2}];
                                end
                            end

                            aux_val = thisOffset.item(0);
                            next_index_base = aux_val.getNodeIndex;

                        end

                    end

                    pb.checkpoint();

%                     % to debug not supported components
%                     if(~bCompParsedOk)
%                         disp_string_framed(2, sprintf('Component %d from %s could not be parsed', comp_idx, filename) );
% 
%                         fprintf(2, sprintf('%s\n', disp_xml_nodes(thisValue) ) );
%                         fprintf(2, sprintf('%s\n', disp_xml_nodes(thisCode) ) );
% 
%                     end

                end

            end
            
            comp_idx = comp_idx + 1;
            
            if(bDesktop)
                pb.end_loop();
            end
            
        end
    
        annSet_idx = annSet_idx + 1;
        
    end
    
    if(bDesktop)            
        % destroy the progress bar.
        pb.delete;
    end
    
    if( bAnnSetFound )
        
        % Assume the heartbeat fiducial point in the dominant wave of the QRS
        % complex
        for ii = 1:length(single_lead_positions)
            if( isempty(single_lead_positions(ii).R) )
                aux_val = [];
            else
                aux_val = single_lead_positions(ii).R(~isnan(single_lead_positions(ii).R));
            end

            if( ~isempty(single_lead_positions(ii).Q) )
                bAux = ~isnan(single_lead_positions(ii).Q) & isnan(aux_val);
                aux_val(bAux) = single_lead_positions(ii).Q(bAux);
            end

            if( ~isempty(single_lead_positions(ii).S) )
                bAux = ~isnan(single_lead_positions(ii).S) & isnan(aux_val);
                aux_val(bAux) = single_lead_positions(ii).S(bAux);
            end

            single_lead_positions(ii).qrs = aux_val;
        end

        % filter repetitions, one label per heartbeat.
        [ann.time, aux_idx] = unique_w_tolerance(ann.time, round(0.15*heasig.freq) );
        ann.anntyp = ann.anntyp(aux_idx);

        % filter repetitions, one wave fiducial point per heartbeat.
        for ii = 1:length(single_lead_positions)
            for fname = rowvec(fieldnames(single_lead_positions(ii)))
                single_lead_positions(ii).(fname{1}) = unique_w_tolerance(single_lead_positions(ii).(fname{1}), round(0.15*heasig.freq) );
            end
        end

        single_lead_positions = groupnlineup_waves(single_lead_positions);

    end
    
end
