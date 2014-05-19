%GETPRIOR Get class prior probabilities of dataset
%
%   [PRIOR,LABLIST] = GETPRIOR(A,WARNING)
%
% INPUT
%   A        Dataset
%   WARNING  1: Generate warning if priors are not set and should be
%               computed from class frequencies (default)
%   WARNING  0: Suppress warning message
%
% OUTPUT
%   PRIOR    Class prior probabilities
%   LABLIST  Label list
%
% DESCRIPTION
% Returns the class prior probabilities as defined in the dataset A.
% In LABLIST the corresponding class labels are returned.
%
% Note that if these are not set (A.PRIOR = []), the class frequencies
% are measured and returned. Use ISEMPTY(A,'prior') to test whether 
% A.PRIOR = [].
%
% If A has soft labels, these are used to estimate the class frequencies. 
% If A has target labels, an error is returned since in that case, no 
% classes are defined.
%
% SEE ALSO
% DATASETS, SETPRIOR
