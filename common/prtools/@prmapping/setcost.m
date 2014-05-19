%SETCOST Reset classification cost matrix of mapping
%
%   W = SETCOST(W,COST,LABLIST)
%
% The classification cost matrix of the dataset W is reset to COST.
% W has to be a trained classifier. COST should have size [C,C+1],
% if C is the number of classes assigned by W.
% COST(I,J) are the costs of classifying an object of class I
% as class J. Column C+1 generates an alternative reject class and
% may be omitted, yielding a size of [C,C].
% An empty cost matrix, COST = [] (default) is interpreted as
% COST = ONES(C) - EYE(C) (identical costs of misclassification).
%
% In LABLIST the corresponding class labels may be supplied.
% LABLIST may have only class names of the existing labels assigned
% by W, stored in W.LABELS.
