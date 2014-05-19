%VERTCAT Vertical concatenation of datafiles (object extension)
%
%    C = [A;B]
%
% The datafiles A and B are vertically concatenated, i.e. the
% objects of B are added to the dataset A. This is consistent with 
% the vertical concatenation of datasets.
%
% The datafiles should be of the same type. It is assumed that
% the preprocessing and postprocessing fields are equal. If not, a
% warning is generated and those of A are used.
%
% See also DATAFILES
