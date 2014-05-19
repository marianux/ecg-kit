%GETIDENT Get fields of object descriptors
%
%    IDENTFIELD = GETIDENT(A,FIELD,L)
%
% INPUT
%   A          Dataset
%   FIELD      Character string: name of structure field of IDENT field of A.
%   L          Vector of indices pointing to desired objects in A, 
%              default: all.
%
% OUTPUT
%   IDENTFIELD Cell array of size (L,1) containing the requested field of the
%              objects A(L,:).
%
% If FIELD is the empty string ('') the entire ident structure is returned.
% If the requested field does not exist IDENTFIELD = [];
%
% Note the ident field of datasets was originally intended for an
% identification of the individual objects. Later its usage was extended 
% to a field for storing general information on objects. For that reason 'old'
% datasets without a structure in the ident field are transformed such that
% this information is stored in a subfield IDENT in the ident field. It can
% be retrieved by GETIDENT(A) or GETIDENT(A,J).
%
% IDENTFIELD is a cell array as arbitrary parameters may be stored. If
% these are doubles, e.g. after A = SETIDENT(A,[1:SIZE(A,1)]'), they can
% be easily converted by N = CELL2MAT(GETIDENT(A));
%
% For backward compatibility the following holds: If FIELD = 'string'
% then IDENTFIELD contains a character array of the object identifiers stored
% in A.IDENT.IDENT. If these are integers they are converted to strings.
