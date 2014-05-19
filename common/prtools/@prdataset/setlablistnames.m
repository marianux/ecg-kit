%SETLABLISTNAMES Reset the names of label lists
%
%       A = SETLABLISTNAMES(A,NEW_NAMES,J)
%       A = SETLABLISTNAMES(A,NEW_NAMES,OLD_NAMES)
%
% INPUT
%   A          - Dataset
%   NEW_NAMES  - String or character array
%   J          - Vector identifying lablist numbers
%   OLD_NAMES  - String or character array
%
% OUTPUT
%   A          - Dataset
%
% DESCRIPTION
% One or more names of the label lists are reset.
% J should have as many elements as names stored in NEW_NAMES.
% Default J=1, in which case NEW_NAMES should be a string.
% In case existing names are identified by OLD_NAMES, this
% character array should have as many names as NEW_NAMES.
%
% SEE ALSO
% DATASETS, MULTI_LABELING, ADDLABELS, ADDLABLIST, CHANGELABLIST,
% CURLABLIST, GETLABLISTNAMES
