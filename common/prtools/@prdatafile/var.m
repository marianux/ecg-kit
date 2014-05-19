%VAR Datafile overload
%
%   [V,U] = VAR(A,W)
%
% A - Datafile
% W - Vector of weights, one per object
%
% V - Vector of feature variances
% U - Mean vector
%
% Computes variance V and mean U in a single run for speed.
% Objects are assumed to have the same number of features.
% Take care that the feature size of A has been correctly set.
% The routine is useful in case the data is too large to be
% converted to a dataset first.
