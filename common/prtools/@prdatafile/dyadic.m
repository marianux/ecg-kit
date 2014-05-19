%DYADIC Dyadic datafile operations
%
%  C = DYADIC(A,P,B,Q)
%
% Computes C = P*A + Q*B
%
% This datafile operation is, like others, either stored as
% a preprocessing or as a postprocessing for datafiles using
% a call to DYADICM.
%
% Note that in P a function name can be stored and in Q a cell
% array with a set of parameters. In that case effectively 
% feval(p,a1,a2,q{:}) will be executed inside DYADICM.
%
% SEE ALSO
% DATAFILES, DYADICM
