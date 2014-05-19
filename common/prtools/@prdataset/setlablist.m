%SETLABLIST Set names of classes or targets
%
%    A = SETLABLIST(A,LABLIST)
%
% LABLIST should be a column vector of C elements (integers, characters,
% or strings), in which C is the number of classes or targets of A.
% In case of multiple label lists this resets the current label list.
%
%    A = SETLABLIST(A)
%
% Removes entries in the lablist of A to which no objects are assigned,
% i.e. it remove empty classes. This command also removes duplicates in the
% lablist. An example of the merge of two classes in a 3-class dataset X
% with class names 'A', 'B' and 'C' can be realized by
%   X = setlablist(X,char('A','B','A'));
%   X = setlablist(X);
%
% SEE ALSO MULTI_LABELING
