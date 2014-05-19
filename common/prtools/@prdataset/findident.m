%FINDIDENT Determine indices of objects having specified identifiers
%
%   J = FINDIDENT(A,IDENT,FIELD)
%
% INPUT
%   A      Dataset
%   IDENT  Object identifiers, see SETIDENT
%   FIELD  Desired field, default 'IDENT'.
%
% If IDENT is a set of object identifiers then J is a vector with indices
% to the first object in A that matches IDENT in the specified FIELD.
% If IDENT is a single object identifier then J is a set of indices to
% all objects in A having the given identifier.
