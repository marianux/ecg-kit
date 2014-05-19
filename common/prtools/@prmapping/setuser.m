%SETUSER Set user field in mapping
%
%   W = SETUSER(W,S,FIELD)
%
% INPUT
%   W      Mapping
%   S      Variable to be stored in the user field 
%   FIELD  Desired field, default 'USER'.
%
% OUTPUT
%   W      Updated mapping
%
% DESCRIPTION
% Set or reset the subfield FIELD of the user field of W by S.
%
% Note the the USER field of mappings was originally intended for a user
% defined description of mappings. Later its usage was extended to a field
% for storing general information on mappings. For that reason 'old'
% mappings without a structure in the user field are transformed such that
% this information is stored in a subfield USER in the user field. It can
% be retrieved by GETUSER(A,'USER').
%
% Note also that for reasons of backward compatibility the parameter order
% of the SETUSER command differs from similar Matlab commands like
% SETFIELD: first field content, then field name.
