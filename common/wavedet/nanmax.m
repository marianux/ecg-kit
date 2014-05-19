function [varargout] = nanmax(varargin)
%NANMAX Maximum value, ignoring NaNs.
%   M = NANMAX(A) returns the maximum of A with NaNs treated as missing. 
%   For vectors, M is the largest non-NaN element in A.  For matrices, M is
%   a row vector containing the maximum non-NaN element from each column.
%   For N-D arrays, NANMAX operates along the first non-singleton
%   dimension.
%
%   [M,NDX] = NANMAX(A) returns the indices of the maximum values in A.  If
%   the values along the first non-singleton dimension contain more than
%   one maximal element, the index of the first one is returned.
%  
%   M = NANMAX(A,B) returns an array the same size as A and B with the
%   largest elements taken from A or B.  Either one can be a scalar.
%
%   [M,NDX] = NANMAX(A,[],DIM) operates along the dimension DIM.
%
%   See also MAX, NANMIN, NANMEAN, NANMEDIAN, NANMIN, NANVAR, NANSTD.

%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 1.1.8.1 $  $Date: 2010/03/16 00:15:49 $

% Call [m,ndx] = max(a,b) with as many inputs and outputs as needed
[varargout{1:nargout}] = max(varargin{:});
