%DELLABLIST Delete a label list from dataset
%
%	B = DELLABLIST(A,LABLISTNAME)
%	B = DELLABLIST(A,LABLISTNUMBER)
%
% INPUT
%   A              - Dataset
%   LABLISTNAME    - String to identify the label list to be deleted
%   LABLISTNUMBER  - Number to identify the label list to be deleted
%
% OUTPUT
%   B              - Dataset
%
% DESCRIPTION
% In the multi-label system for datasets, additional labellings can be
% added by ADDLABELS. This is stored in the LABLIST and NLAB fields of
% the dataset. By this command, DELLABLIST such a labeling can be removed.
%
% SEE ALSO
% DATASETS, MULTI_LABELING, ADDLABELS, ADDLABLIST, CHANGELABLIST, CURLABLIST
