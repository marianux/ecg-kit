%GETSIZE Dataset size and number of classes
%
%  [M,K,C] = GETSIZE(A,DIM)
%
% INPUT
%   A    Dataset
%   DIM  1,2 or 3 : the number of the output argument to be returned
%
% OUTPUT
%   M    Number of objects
%   K    Number of features
%   C    Number of classes
%
% DESCRIPTION
% Returns size of the dataset A and the number of classes. C is determined
% from the number of labels stored in A.LABLIST. If DIM = 1,2 or 3, just 
% one of these numbers is returned, e.g. C = GETSIZE(A,3).
