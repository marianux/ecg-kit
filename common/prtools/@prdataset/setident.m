%SETIDENT Set object identifiers
%
%   A = SETIDENT(A,IDENT,FIELD,L)
%
% INPUT
%   A      Dataset
%   IDENT  Object identifiers, size (N,K)
%   FIELD  Desired field, default 'IDENT'.
%   L      Vector of indices of objects to be updated (optional; default: all)
%          length(L) = N.
%
% OUTPUT
%   A      Updated dataset
%
% DESCRIPTION
% Set or reset the subfield FIELD of the ident field of A by IDENT.
% IDENT should be an array of size (N,K), with arbitrary K. 
%
% Note the ident field of datasets was originally intended for an
% identification of the individual objects. Later its usage was extended 
% to a field for storing general information on objects. For that reason 'old'
% datasets without a structure in the ident field are transformed such that
% this information is stored in a subfield IDENT in the ident field. It can
% be retrieved by GETIDENT(A,'IDENT').
%
% The default FIELD is 'IDENT'. To reset the entire IDENT give
% A = SETIDENT(A,IDENT,''), in which IDENT is a structure array of the
% right size and including a subfield named also IDENT.
%
% The new structure is checked or created by A = SETIDENT(A);
%
% Note also that for reasons of backward compatibility the parameter order
% of the SETIDENT command differs from similar Matlab commands like
% SETFIELD: first field content, then field name.
