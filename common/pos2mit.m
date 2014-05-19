function pos2mit(ann_filename, ann_samples, ann_type)
% output = pos2mit(ann_samples,ann_label,ann_lead,ann_path,ann_filename)
% Function to create a MIT annotator file, with any anots and they are
% considered as "normal" beats with the aim of being usable in bsb
% interface mode.
%
% Input:
%    ann_samples: marks of annots to convert.
%    ann_label: string descripting the marks (i.e. Toff) (by default:
%    'anot1')
%    ann_lead: string descripting from which lead these marks belong. (by
%    default: '1')
%    ann_path: path where annotation file will be saved. (by default: pwd)
%    ann_filename: name of the signal recording from annotations come out. (by
%    default: 'anot2mit')
%
% NO UPDATES

ann_samples = colvec(ann_samples);

ann_samples = ann_samples(~isnan(ann_samples));
s.time = ann_samples;

if( nargin < 3 || isempty(ann_type))
    s.anntyp(1:length(s.time),1) = 'N';
else
    s.anntyp = ann_type;
end

s.subtyp(1:length(s.time),1) = ' ';
s.chan = s.subtyp;
s.num = s.chan;
s.aux(1:length(s.time),1) = ' ';
writeannot(ann_filename,s);
