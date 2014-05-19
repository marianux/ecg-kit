%GETDATA Retrieve data field from mapping or dataset
%
%    DATA = GETDATA(W)
%    DATA = GETDATA(W,N)
%    DATA = GETDATA(W,FIELD)
%
%    [DATA,NLAB,IDENT] = GETDATA(A,CLASSES)
%
% INPUT
%   W        Mapping
%   N        Index of cell array to be used for DATA (optional)
%   FIELD    String, name of data field structure to be retrieved
%            (optional)
%
%   A        Dataset
%   CLASSES  Indices of classes to be retrieved (optional)
%
% OUTPUT
%   DATA     Retrieved data
%   NLAB     Numeric labels of retrieved objects
%   IDENT    Object identifiers of retrieved objects 
%
% DESCRIPTION
% This routine can be used to retrieve data stored in the data fields of a
% mapping W or a dataset A.
%
% The data field of a mapping W consists of a matrix, cell array or a
% structure. If neither N nor FIELD are given, the entire data field is
% returned in DATA. If the data field is a structure, the string FIELD may
% point to the desired structure field to be returned. In case it is a cell
% array N may point to the desired element.
%
% The data field of a dataset A is an array of doubles. If the index vector 
% CLASSES is given just the data of the objects (rows of A.DATA) that belong 
% to the corresponding classes are returned. CLASSES are the class numbers
% of the classes defined in the current LABLIST. 
%
% By default all data is returned. In NLAB the class numbers of the
% selected objects are given, these are the indices to the label list that
% can be retrieved by LABLIST = GETLABLIST(A).
%
% IDENT returns the vector or cell array with the corresponding object
% identifiers stored in the IDENT field os the dataset A.
%
% SEE ALSO
% MAPPINGS, DATASETS
