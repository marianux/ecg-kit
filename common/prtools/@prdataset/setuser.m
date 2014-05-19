%SETUSER Set the user field of a dataset
%
%   A = SETUSER(A,S,FIELD)
%
% INPUT
%   A      Dataset
%   S      Variable to be stored in the user field 
%   FIELD  Desired field, default 'USER'.
%
% OUTPUT
%   A      Updated dataset
%
% DESCRIPTION
% Set or reset the subfield FIELD of the user field of A by S.
%
% Note the the USER field of datasets was originally intended for a user
% defined description of datasets. Later its usage was extended to a field
% for storing general information on datasets. For that reason 'old'
% datasets without a structure in the user field are transformed such that
% this information is stored in a subfield USER in the user field. It can
% be retrieved by GETUSER(A,'USER').
%
% Note also that for reasons of backward compatibility the parameter order
% of the SETUSER command differs from similar Matlab commands like
% SETFIELD: first field content, then field name.
