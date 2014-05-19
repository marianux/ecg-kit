%MULTI_LABELING Info on the PRTools multiple labeling system.
%
% This is not a command, just an information file.
%
% In the PRTools datasets the data is stored together with the class labels
% for the objects. For the crisp and soft labels, as well as for targets a 
% system of multiple labels has been constructed by which several label sets, 
% the corresponding class names, prior probabilities and classification costs
% can be stored. There is always one of them active. 
%
% This system, available from PRTools4.1 is build on top of the old, single 
% labling system used until PRTools4.0. It is implemented in such a way
% that the existing software is not be disturbed by it. The following new
% commands are available for using multiple labeling. Recall that PRTools
% datasets do not store the object labels as such, but as a label list
% (list of class names), stored in the dataset LABLIST field and by a vector
% of LABLIST indices (one index per object) stored in the NLAB field.
%
% - [B,N] = ADDLABELS(A,LABELS,LABLISTNAME) :
% This adds an additional set of labels to the dataset A, stores it as a
% label list using the given label list name (LABLISTNAME). This name, or the
% returned number N may be used to activate the labels in future. It is
% also activated now.
%
% - B = CHANGELABLIST(A,LABLISTNAME) :
%   B = CHANGELABLIST(A,N)
% Activates a previously stored label list.
%
% - B = DELLABLIST(A,LABLISTNAME) :
%   B = DELLABLIST(A,N) :
% Deletes an existing label list.
%
% - GETLABLISTNAMES(A) :
% Returns the names of the stored label lists.
%
% - SETLABLISTNAMES(A,NEW_NAMES,OLD_NAMES) :
% Change the names of existing label lists.
%
% - [N,LABLISTNAME] = CURLABLIST(A) :
% Returns the number and names of the current label list.
%
% - [B,N] = ADDLABLIST(A,LABLIST,LABLISTNAME) :
% This is a low level routine, of interest of PRTools maintenance
% programmers only. It organises the storage of information. This routine
% just stores a new list of label names and makes space for the NLAB
% indices. It should thereby by followed by a SETNLAB call. Read the
% ADDLABLIST help file for more information on the implementation.
%
% The NLAB dataset field is extended to a MxN array, in which M is the
% number of objects and N is the number of label lists. The existing
% routines like GETLABELS, SETLABELS, GETNLAB, SETNLAB, GETLABLIST,
% SETLABLIST, GETPRIOR, SETPRIOR, GETCOST and SETCOST refer to the current
% label list only. Interaction with the multiple label system is done
% automatically.
%
% SEE ALSO
% ADDLABELS, CHANGELABLIST, DELLABLIST, GETLABLISTNAMES, SETLABLISTNAMES,
% ADDLABLIST, GETLABELS, SETLABELS, GETNLAB, SETNLAB, GETLABLIST,
% SETLABLIST, GETPRIOR, SETPRIOR, GETCOST, SETCOST

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

