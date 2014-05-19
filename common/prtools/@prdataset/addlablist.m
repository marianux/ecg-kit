%ADDLABLIST Add another label list to dataset definition
%
% [B,N] = ADDLABLIST(A,LABLIST,LABLISTNAME)
%
% INPUT
%   A            - Dataset
%   LABLIST      - Label list to be added
%   LABLISTNAME  - Optional name for this label list
%
% OUTPUT
%   B            - Dataset
%   N            - Number for the new label list
%
% DESCRIPTION
% This command adds an additional label list to a dataset. It also
% activates this label list as the current one. LABLISTNAME may be
% given to identify this label list in the CHANGELABLIST command.
% Alternatively, the label list number may be used there.
%
% Multiple label lists are only supported for crisp label types and 
% not for the soft and target label types. See SETLABELS.
%
% In addition to storing the new label list the following changes in
% the dataset are made:
% - If the dataset did not contain a mulitple label setup yet, it is
%   implemented. The following definition is used for L label lists:
%   - B.LABLIST becomes a cell array of size (L+1,2)
%   - B.LABLIST(I,1) contains label list I
%   - B.LABLIST(I,2) contains the corresponding prior probabilities
%     (to be substituted in B.PRIOR by CHANGELABLIST)
%   - B.LABLIST(I,3) contains the corresponding error costs
%     (to be substituted in B.COST by CHANGELABLIST)
%   - B.LABLIST(I,4) contains the corresponding label type (crisp, soft or
%     targets).
%   - B.LABLIST(L+1,1) contains a char array with the label list names
%     These names can be retrieved by GETLABLISTNAMES
%   - B.LABLIST(L+1,2) contains the index of the current label list.
%   - B.LABLIST(L+1,3) contains a vector with target sizes (number of
%   	columns) for each label list. Used by SETTARGETS and GETTARGETS.
%   - B.NLAB is an array of size (M,L), in which M is the number of
%     objects in B
% - The number of columns in NLAB is extended by one. This columns is
%   filled by zeros (no labels defined yet).
% - B.PRIOR is made empty
% - B.COST is made empty
% - An new cell pair is added in B.LABLIST just before the last pair.
%   it is filled by the new label list and an empty prior definition
% - The set of label lists in B.LABLIST{L+1,1) is updated with the new
%   label list name. 
% - The current label list in B.LABLIST{L+1,2) is set to L.
%
% The ADDLABLIST command should be followed by a SETNLAB to store the
% proper new labels for the objects and by a SETPRIOR to update the
% prior probability setting. These actions may be combined by the
% ADDLABELS command. An existing label list may be removed by the
% DELLABLIST command.
%
% [A,N,T0,T1] = ADDLABLIST(A)
%
% In this case just the current label list is returned in N and the start
% and endpoint of the corresponding columns in A.TARGET. If T0<T1: no
% targets are set.
%
% SEE ALSO
% DATASETS, MULTI_LABELING, ADDLABELS, CHANGELABLIST, CURLABLIST, DELLABLIST
