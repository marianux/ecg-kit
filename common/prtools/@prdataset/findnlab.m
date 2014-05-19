%FINDNLAB Determine indices of specified classes (numeric)
%
%    J = FINDNLAB(A,NLAB)
%
% INPUT
%     A        Dataset
%     NLAB     vector with numerical indices of reqested classes
%
% OUTPUT
%     J        vector of sample indices
%
% DESCRIPTION
% J is a column vector of indices to the objects in A that have one of the
% numeric labels stored in NLAB. Note that NLAB should contain indices to
% the labellist of A, which may be retrieved by LABLIST = GETLABLIST(A).
% This command thereby retrieves the indices to all objects a specific class
% or set of classes.
