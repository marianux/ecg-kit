%CHANGELABLIST Change current label list 
%
% B = CHANGELABLIST(A,LABLISTNAME)
% B = CHANGELABLIST(A,LABLISTNUMBER)
%
% INPUT
%   A             - Dataset
%   LABLISTNAME   - Name of desired label list
%   LABLISTNUMBER - Number of desired label list (default 1)
%
% OUTPUT
%   B             - Dataset
%
% DESCRIPTION
% This command makes another label list, stored already in A by
% ADDLABLIST, the current one.
% The default label list is the first one set for A.
%
% SEE ALSO
% DATASETS, MULTI_LABELING, ADDLABLIST, CURLABLIST
