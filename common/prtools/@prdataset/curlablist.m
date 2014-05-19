%CURLABLIST Get current label list 
%
% [LABLISTNUMBER,LABLISTNAME,T0,T1] = CURLABLIST(A)
%
% INPUT
%   A             - Dataset
%
% OUTPUT
%   LABLISTNUMBER - Number of current lablist
%   LABLISTNAME   - Name of current label list
%		T0            - Start column current targets
%   T1            - End column current targets
%
% DESCRIPTION
% The number and name of the current label list are returned.
% If T1<T0: no tergets (or soft labels) set for current label list.
%
% SEE ALSO
% DATASETS, MULTI_LABELING, ADDLABLIST, CHANGELABLIST
