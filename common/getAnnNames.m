%% (Internal) Get names of annotations from annotation structure
%   
%   [AnnNames, all_annotations] = getAnnNames(aux_struct)
% 
% Arguments:
% 
%      + aux_struct: 
% 
%      + retries: times to check the existence
%             
% Output:
% 
%      + AnnNames:
% 
%      + all_annotations:
% 
% Example:
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
function [AnnNames, all_annotations] = getAnnNames(aux_struct)

    AnnNames = [];

    for fname = rowvec(fieldnames(aux_struct))
        if( isfield(aux_struct.(fname{1}), 'time') )
            AnnNames = [AnnNames; cellstr(fname{1}) cellstr('time')];
        end
        if( isfield(aux_struct.(fname{1}), 'qrs') )
            AnnNames = [AnnNames; cellstr(fname{1}) cellstr('qrs')];
        end
    end

    cant_anns = size(AnnNames,1);

    all_annotations = cell(cant_anns,1);
    for ii = 1:cant_anns
        all_annotations{ii} = aux_struct.(AnnNames{ii,1}).(AnnNames{ii,2});
    end
    