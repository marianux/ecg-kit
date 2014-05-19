%SETFEATSIZE Set the feature size of a dataset
%
%   A = SETFEATSIZE(A,FEATSIZE)
%
% INPUT 
%   A         Dataset
%   FEATSIZE  Feature size vector
% 
% OUTPUT
%   A         Dataset
%
% DESCRIPTION
% By default the feature size of a dataset is its number of features, i.e.
% the number of columns in the DATA field of A. If the features are samples
% of a multi-dimensional data item, e.g. the pixels of an image, the
% original size of this data item may be stored in FEATSIZE. The product of
% all elements in FEATSIZE has to be equal to the number of columns in the
% DATA field of A.
