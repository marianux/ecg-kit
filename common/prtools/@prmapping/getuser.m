%GETUSER Return the user field of a mapping
%
%    USERFIELD = GETUSER(W,FIELD)
%
% INPUT
%   W          Mapping
%   FIELD      Character string: name of structure field of USER field of W.
%
% OUTPUT
%   USERFIELD  Requested field
%
% If FIELD is is the empty string ('') the entire user field is returned. 
%
% If the requested field does not exist USERFIELD = [];
%
% Note the the USER field of mappings was originally intended for a user
% defined description of mappings. Later its usage was extended to a field
% for storing general information on mappings. For that reason 'old'
% mappings without a structure in the user field are transformed such that
% this information is stored in a subfield USER in the user field. It can
% be retrieved by GETUSER(W).
