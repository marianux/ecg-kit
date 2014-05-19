%SETNLAB Set label indices (numeric labels) directly
%
%    A = SETNLAB(A,NLAB,J)
%
% Resets the values of A.NLAB(J,N) to NLAB, in which N points
% to the current label list for A. Note that NLAB gives indices
% to the labels defined in this label list. Values of NLAB <= 0 are
% stored but treated as missing labels. For crisp labels only.
%
% SEE ALSO MULTI_LABELING
