%GETCLASSI Get class index
%
%   INDEX = GETCLASSI(A,LABEL)
%
% INPUT
%   A     Dataset
%   LABEL Label used to label the class in A
%
% OUTPUT
%   INDEX Index of LABEL in the label list of A
%
% DESCRIPTION
% In some routines like SELDAT and ROC classes should be defined
% by their index in the label list of the dataset. This label list
% can be retrieved by GETLABLIST or CLASSNAMES. GETCLASSI can be used to
% find the index of the label directly.
%
% SEE ALSO
% DATASETS, SETLABELS, SETLABLIST GETLABLIST, CLASSNAMES
