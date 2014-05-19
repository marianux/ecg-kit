%ADDLABELS Add a new labeling to an existing dataset
%
%	[B,LABLISTNUMBER] = ADDLABELS(A,LABELS,LABLISTNAME)
%
% INPUT
%   A              - Dataset
%   LABELS         - Vector or character array with labels
%   LABLISTNAME    - String to identify the new labeling
%
% OUTPUT
%   B              - Dataset
%   LABLISTNUMBER  - Number for the new label list
%
% DESCRIPTION
% This adds a new labeling to the given dataset A. It has thereby a
% multiple labeling. The new labeling is immediately activated and made the
% current one. See MULTI_LABELING for a description of the multiple labeling
% system. See ADDLABLIST for implementation details.
% Desired labellings may be set by CHANGELABLIST using the LABLISTNUMBER 
% or LABLISTNAME. The original, first labeling of a dataset has 
% LABLISTNUMBER = 1 and LABLISTNAME = 'default'.
%
% Use SETPRIOR and SETCOST to set priors and costs for the new labeling.
% A defined labeling can be removed by DELLABLIST.
% This multiple labeling system is implemented for crisp labels only.
%
% SEE ALSO
% DATASETS, MULTI_LABELING, ADDLABLIST, CHANGELABLIST, DELLABLIST, SETPRIOR
