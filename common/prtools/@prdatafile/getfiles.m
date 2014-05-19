%GETFILES Get datafile FILES field
%
%	   FILES = GETFILES(A)
%    FILE  = GETFILES(A,N)
%
% Retrieves the filenames on which the datafile A is defined,
% or just the name of the N-th file.
% If N larger than the number of stored files, FILE = [].
%
% Note that if A.RAW is true some or all FILES may be directories.
% If A.RAW is false, all FILES are mat-files containing datasets.
