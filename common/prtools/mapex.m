%MAPEX Train and execute arbitrary untrained mapping
%
%    B = A*(W*MAPEX) = A*(A*W)
%    B = A*MAPEX(W,PAR1,PAR2, ...) 
%    B = A*MAPEX(UMAP,PAR1,PAR2, ...)
%
% INPUT
%    A        Dataset or datafile or double
%    W        Untrained mapping
%    UMAP     String with name of untrained mapping
%    PAR1     Parameters of untrained mapping W or UMAP
%
% OUTPUT
%    B               Resulting dataset, datafile or double array
%
% DESCRIPTION
% This routine facilitates the construction of shortcuts by training and
% executing an untrained mapping by the same data. It is typically useful
% for mappings like PROXM and SCALEM that are often trained and exectuted
% with the same data.
%
% EXAMPLE
% mink1 = proxm('m',1)*mapex
% mink1 = mapex(proxm,'m',1)    % the same
% mink1 = mapex('proxm','m',1)  % the same
% Herewith D = A*mink1 computes a dissimilarity matrix based on the
% Minkowsky_1 metric between all objects in A.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, PROXM, SCALEM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function out = mapex(p1,varargin)

    
  if nargin == 0
    out = prmapping(mfilename,'combiner');
  elseif nargin == 1 && ismapping(p1)
    out = prmapping(mfilename,'fixed',p1);
  elseif nargin >= 1 && ischar(p1)
    out = prmapping(mfilename,'fixed',{p1,varargin{:}});
  elseif nargin > 1 && ismapping(p1)
    out = prmapping(mfilename,'fixed',{getmapping_file(p1),varargin{:}});
  elseif (isdataset(p1) || isa(p1,'double')) && ismapping(varargin{1})
    out = p1*(p1*varargin{1});
  elseif (isdataset(p1) || isa(p1,'double')) && nargin == 2
    out = p1*feval(varargin{1},p1);
  elseif (isdataset(p1) || isa(p1,'double')) && nargin > 1
    out = p1*feval(varargin{1},p1,varargin{2:end});
  else
    error('Illegal input')
  end