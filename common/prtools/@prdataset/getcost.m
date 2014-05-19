%GETCOST Get classification cost matrix
%
%  [COST,LABLIST] = GETCOST(A)
%
% Returns the classification cost matrix as defined for the dataset A.
% An empty cost matrix is interpreted as equal costs for misclassification.
% In that case COST = ONES(C+1) - EYE(C+1) is returned, if C is the number
% of classes. Row C+1 and column C+1 in COST refer to unlabeld objects.
% In LABLIST the class labels are returned.
%
% If A has target labels an error is returned as in that case no classes
% and thereby no costs are defined.
%
% See DATASETS
