%GETTARGETS Get targets of dataset
%
%    [TARGETS,LABLIST] = GETTARGETS(A)
%
% Gets the targets (or soft labels) of the objects in the dataset A.  If A
% has crisp labels they are converted to targets first. See SETLABTYPE for
% conversion rules.
%
%    [TARGETS,LABLIST] = GETTARGETS(A,'soft')
%
% Forces the return of soft labels. This command is identical to
% GETLABELS(A,'soft'). If A has crisp labels or target labels they are
% converted to soft labels first. See SETLABTYPE for conversion rules.
%
% See also SETLABTYPE, GETLABELS, SETLABELS, SETTARGETS
