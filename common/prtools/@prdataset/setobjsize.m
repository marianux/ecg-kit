%SETOBJSIZE Reset the object size of a dataset
%
%    A = SETOBJSIZE(A,OBJSIZE)
%
% By default, the object size of a dataset A is given by the number of objects, 
% i.e. the number of rows in the DATA field of A. If the objects are samples of
% a multi-dimensional data item, e.g. the pixels of an image, the original size
% of this data item may be stored in OBJSIZE. The product of all elements in 
% OBJSIZE has to be equal to the number of rows in the DATA field of A.
