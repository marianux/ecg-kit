%STD Datafile overload
%
%   [S,U] = STD(A,FLAG,DIM)
%
% A - Datafile
% FLAG - FLAG==0 for the default normalization by N-1, or
%        FLAG==1 for normalizing by N.
%
% S - Vector of feature standard deviatons
% U - Mean vector
%
% Computes std. dev. S and mean U in a single run for speed.
% Objects are assumed to have the same number of features.
% Take care that the feature size of A has been correctly set.
% The routine is useful in case the data is too large to be
% converted to a dataset first.
