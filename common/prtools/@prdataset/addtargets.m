%ADDTARGETS Add a new targets to an existing dataset
%
%	[B,LABLISTNUMBER] = ADDTARGETS(A,TARGETS,LABLISTNAME)
%
% INPUT
%   A              - Dataset
%   TARGETS        - Array or dataset with targets
%   LABLISTNAME    - String to identify the new set of targets
%
% OUTPUT
%   B              - Dataset
%   LABLISTNUMBER  - Number for the new label list
%
% DESCRIPTION
% This adds a new set of targets to the given dataset A. It has thereby a
% multiple labeling and/or targets. The new targets are immediately activated
% and made the current one. See MULTI_LABELING for a description of the 
% multiple labeling/target system. See ADDLABLIST for implementation details.
% Desired target sets or labellings may be set by CHANGELABLIST using the 
% LABLISTNUMBER or LABLISTNAME. The original, first labeling of a dataset has 
% LABLISTNUMBER = 1 and LABLISTNAME = 'default'.
%
% A defined target set can be removed by DELLABLIST.
%
% SEE ALSO
% DATASETS, MULTI_LABELING, ADDLABLIST, ADDLABELS, CHANGELABLIST, DELLABLIST
