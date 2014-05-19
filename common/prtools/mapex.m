%MAPEX Train and execute arbitrary untrained mapping
%
%    B = A*(W*MAPEX) = A*(A*W)
%
% INPUT
%    A        Dataset or datafile or double
%    W        Untrained mapping
%
% OUTPUT
%    B               Resulting dataset, datafile or double array
%
% DESCRIPTION
% This routine facilitates the construction of shortcuts by training and
% executing a mapping by the same data.
%
% EXAMPLE
% mink1 = proxm([],'m',1)*mapex
% Herewith D = A*mink1 computes a dissimilarity matrix based on the
% Minkowsky_1 metric between all objects in A.
%
% DATASETS, MAPPINGS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function out = mapex(p1,p2)

    
  if nargin == 0
    out = prmapping(mfilename,'combiner');
  elseif ismapping(p1)
    out = prmapping(mfilename,'trained',p1);
  elseif isdataset(p1) && ismapping(p2)
    w = getdata(p2);
    out = p1*(p1*w);
  else
    error('Illegal input')
  end