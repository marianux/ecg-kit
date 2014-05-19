%SETCOST Reset classification cost matrix of dataset
%
%   A = SETCOST(A,COST,LABLIST)
%
% The classification cost matrix of the dataset A is reset to COST.
% COST should have size [C,C+n], n >= 0, if C is the number of classes.
% COST(I,J) are the costs of classifying an object of class I
% as class J. Columns C+j generate an alternative reject classes and
% may be omitted, yielding a size of [C,C].
% An empty cost matrix, COST = [] (default) is interpreted as
% COST = ONES(C) - EYE(C) (identical costs of misclassification).
%
% In LABLIST the corresponding class labels may be supplied.
% LABLIST may have only class names of the existing classes in A.
% Reset class names first by SETLABLIST if necessary.
%
% Alternatively, for classification matrices, LABLIST may refer to
% the class names stored in the feature labels. This should be used
% with care as it may disturb the existing lableling of A.
%
% If LABLIST is not given, the order defined by the existing LABLIST
% for A (determined by [NLAB,LABLIST] = renumlab(LABELS)) is used.
