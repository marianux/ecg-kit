%SETFEATDOM Reset feature domains of dataset  
%
%   A = SETFEATDOM(A,FEATDOM,K)
%
% FEATDOM has to be a cell-array, defining for each feature
% the desired domain. The following formats are supported:
% []         : empty, all numbers are allowed
% R size(N,2): Defines N ranges on the real axis, 
%              lowerbounds are R(:,1), upperbounds R(:,2).
%              N SHOULD BE AT LEAST 2. E.g:
%              R = [0 1; 0 1] defines a single region between 0 and 1.
%              R = [-1 1; 5 10; 11 11] defines two regions and
%              a single value (11).
% J size(1,N): Set of N integer values
% S size(N,F): String array of N strings. The feature values in the
%              DATA field of A will be coded by integers from 1 to N.
%              A feature like this may be assigned to the dataset
%              by A(:,j) = S, in which S is a string array, containing
%              a string for each object.  
% By this routine and also by all changes in the dataset, TESTFEATDOM 
% will be called to check consistency between data and feature domain
% descriptions.
% 
% If given, K has to be a set of indices with length equal to
% length(FEATDOM) pointing to the feature domains to be reset.
% Checking is done as well. So reset data first if needed.
