%GETLABELS Get labels or soft labels of a dataset
%
%	  [LABELS,LABLIST] = GETLABELS(A,TYPE,LABLISTNAME)
%
% INPUT
%  A            Dataset
%  TYPE         Label type for conversion, e.g. 'crisp' or 'soft'. 
%               Default: no conversion.
%  LABLISTNAME  Desired lablist, default: present on of A.
%
% OUTPUT
%  LABELS
%  LABLIST
%
% DESCRIPTION
% Get the labels (crisp or soft) of the objects in the dataset A.
% If A has target labels they are converted to soft labels first. 
% See SETLABTYPE for conversion rules.
% LABLIST is the unique set of labels of A and is thereby identical to
% the class names of A.
%
% TYPE = 'soft' forces the return of soft labels after conversion (if 
% necessary). This is identical to GETTARGETS(A,'soft').
% Note that soft labels are not names, but memberships to all classes.
% If A has crisp labels or target labels they are converted to soft
% labels first. See SETLABTYPE for conversion rules.
% 
% SEE ALSO
% SETLABTYPE, GETTARGETS, SETLABELS, SETTARGETS, MULTI_LABELING
