%SETDATA Set data field in mapping or dataset
%
% 	W = SETDATA(W,DATA,FIELD)
% 	W = SETDATA(W,DATA,N)
%
%   A = SETDATA(A,DATA,FEATLAB)
%
% INPUT
%   W        Mapping
%   DATA     Data to be put in the data field.
%   FIELD    String, name of data field structure to be used for DATA
%            (optional)
%   N        Index of cell array to be used for DATA (optional)
%
%   A        Dataset
%   FEATLAB  Desired feature labels for dataset (optional)
%
% OUTPUT
%   W        Mapping
%   A        Dataset
%
% DESCRIPTION
% This routine can be used to store data in the data fields of either a
% mapping W or a dataset A.
%
% The data field of a mapping consists of a matrix, cell array or a structure. 
% In case nor FIELD neither N are given the entire dats field is replaced by 
% DATA. In case the data field is a structure DATA is assigned to field FIELD.
% In case the datafield is a cell array, DATA is assigned to cell N.
% DATA cannot be a structure in case of untrained or fixed mappings.
%
% The data field of a dataset A is an array of doubles. It is replaced by 
% DATA which may be doubles or another dataset, which will be converted to 
% doubles. The numbers of objects in A and DATA should be equal. The feature 
% labels of A are replaced by FEATLAB, or, if not supplied and DATA is a 
% dataset by the feature labels of DATA. Labels and class probabilities of 
% A are preserved.
%
% SEE ALSO
% MAPPINGS, DATASETS
