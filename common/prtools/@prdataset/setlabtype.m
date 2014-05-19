%SETLABTYPE Reset label type of dataset
%
%    A = SETLABTYPE(A,TYPE,LABELS)
%
% The label type of the dataset A is converted to TYPE ('crisp','soft' or
% 'targets'). A conversion of the dataset fields 'nlab', 'lablist' and
% 'targets' is made where necessary. If given, LABELS replaces the labels
% or targets of A.
%
% EXAMPLE
% a = prdataset(rand(10,5)); % create dataset of 10 objects and 5 features
% a = setlabtype(a,'soft',rand(10,1)); % give it random soft labels
