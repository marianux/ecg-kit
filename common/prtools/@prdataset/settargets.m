%SETTARGETS Reset targets or soft labels of dataset
%
%     A = SETTARGETS(A,TARGETS,LABLIST)
%
% The targets or soft label values of the dataset A are reset to TARGETS. 
% This has to be an array or a dataset of size [M,C], in which M is
% the number of objects in A and C is the number of targets or classes.
% The names of these targets or classes can be supplied in the column
% vector LABLIST (numbers, characters or strings) of length C.
% Alternatively, TARGETS may be a dataset with feature labels
% (FEATLAB) equal to LABLIST.
