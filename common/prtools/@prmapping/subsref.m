%SUBSREF Subscript reference overload of mapping
%
% This routine enables constructions like DATA = W.DATA,
% which is similar to DATA = GETDATA(W).
%
% In addition V = W(I,J) is supported for affine transformations.
% It is again an affine mapping using the [I,J] block of the
% rotation matrix and the elements J of the support vector.
%
% For arbitrary mappings just V = W(:,J) is defined by output
% selection: A*V returns just the features J of A*W.
