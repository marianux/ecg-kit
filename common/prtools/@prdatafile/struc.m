%STRUC Datafield structure overload
%
%	[E,D] = STRUC(A)
%
% INPUT
%   A   Datafile
%
% OUTPUT
%   E   Structure conversion of datafile
%   D   Structure conversion of dataset filed of datafile
%
% DESCRIPTION
% As the datafile object is a child of the dataset object the important
% information stored in the dataset field is not directly visible and
% accessible. This commands facilitates this, especially important to
% visualize the two structures by:
%
% STRUC(A)
% 
% without output arguments: it lists the two structures.
