%SETFEATLAB Reset feature labels of a dataset
%
%   A = SETFEATLAB(A,FEATLAB,K)
%
% INPUT
%   A        Dataset
%   FEATLAB  Feature labels
%   K        Vector of indices of feature labels to be reset  
%
% OUTPUT
%   A        Updated dataset
%
% DESCRIPTION
% Set or reset the feature labels of A. The feature labels FEATLAB of the 
% objects in A should be given as a column array of numbers (vector), or
% characters, or as a string array having as many rows as A has features. 
% If given, LENGTH(K) should be equal to LENGTH(FEATLAB).
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, PRDATASET, GETFEATLAB
