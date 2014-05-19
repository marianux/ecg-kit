%SETLABELS Reset labels of dataset or mapping
%
%   A = SETLABELS(A,LABELS,J)
%   W = SETLABELS(W,LABELS,J)
%
% INPUT
%   A          Input dataset, size [M,K]
%   W          Input mapping (classifier), size [K,C]
%   LABELS     Desired labels. M or length(J) labels for a dataset.
%              C or length(J) labels for a mapping.
% OUTPUT
%   A          Dataset with reset labels
%   W          Mapping with reset labels
%
% DESCRIPTION
% The labels of the dataset A are reset by LABELS. If supplied, the
% index vector J defines the objects for wich LABELS applies. If in
% LABELS just a single label is given all the objects defined by J
% are given that label. If LABELS is empty ([]) or NaN all the objects
% defined by J are marked as unlabeled.
%
% If A has soft labels (label type is 'soft') or has no labels but
% targets (label type is 'targets'), these soft labels or targets are
% replaced by LABELS, provided it has the right size.
%
% For soft labels and targets supplied to relabel a dataset, LABELS may be 
% supplied as a dataset of which the data are used for the soft labels or 
% targets and the feature labels are used to set LABLIST of A.
% 
% The labels stored in a classifier mapping W are assigned as feature
% labels of the resulting dataset D in case a dataset B is applied to W:
% D = A*W.
%
% SEE ALSO 
% DATASETS, MAPPINGS, MULTI_LABELING
