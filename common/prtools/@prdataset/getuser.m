%GETUSER Get the user field of a dataset
%
%    USERFIELD = GETUSER(A,FIELD)
%
% INPUT
%   A          Dataset
%   FIELD      Character string: name of structure field of USER field of A.
%
% OUTPUT
%   USERFIELD  Requested field
%
% If FIELD is is the empty string ('') the entire user field is returned. 
%
% If the requested field does not exist USERFIELD = [];
%
% Note the the USER field of datasets was originally intended for a user
% defined description of datasets. Later its usage was extended to a field
% for storing general information on datasets. For that reason 'old'
% datasets without a structure in the user field are transformed such that
% this information is stored in a subfield USER in the user field. It can
% be retrieved by GETUSER(A).
