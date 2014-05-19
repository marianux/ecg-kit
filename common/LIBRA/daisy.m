function result = daisy(x,vtype,metric)

%DAISY returns a matrix containing all the pairwise dissimilarities
% (distances) between observations in the dataset.  The original
% variables may be of mixed types.
%
% The calculation of dissimilarities is explained in:
%   Kaufman, L. and Rousseeuw, P.J. (1990),
%   "Finding groups in data: An introduction to cluster analysis",
%   Wiley-Interscience: New York (Series in Applied Probability and
%   Statistics), ISBN 0-471-87876-6.
%
% Required input arguments:
%       x : Data matrix (rows = observations, columns = variables)
%   vtype : Variable type vector (length equals number of variables)
%           Possible values are 1  Asymmetric binary variable (0/1)
%                               2  Nominal variable (includes symmetric binary)
%                               3  Ordinal variable
%                               4  Interval variable
%  
% Optional input arguments:
%   metric : Metric to be used 
%            Possible values are 'eucli' Euclidian (all interval variables, default)
%                           'manha' Manhattan
%                           'mixed' Mixed (not all interval variables, default)
%
% I/O:
%   result=daisy(x,vtype,'eucli')
%
% Example (subtracted from the referenced book)
%   load flower.mat
%   result=daisy(flower,[2 2 1 2 3 3 4]);
%
% The output of DAISY is a structure containing:
%   result.disv       : dissimilarities (read row by row from the
%                       lower dissimilarity matrix)
%   result.metric     : metric used
%   result.number     : number of observations
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at:
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Guy Brys and Wai Yan Kong (May 2006)

%Checking and filling in the inputs
if (nargin<2)
    error('Two input arguments required')
elseif (nargin<3)
    if (sum(vtype)~=4*size(x,2))
        metric = 'mixed';
        metr = 0;
    else
        metric = 'eucli';
        metr = 1;
    end
elseif (nargin==3)
    if strcmp(metric,'eucli')
        metr=1;
    elseif strcmp(metric,'manha')
        metr=2;
    elseif strcmp(metric,'mixed')
        metr=0;
    end
end

%Standardizing in case of mixed metric
if (sum(vtype)~=4*size(x,2))
    colmin = min(x);
    colextr = max(x)-colmin;
    x = (x - repmat(colmin,size(x,1),1))./repmat(colextr,size(x,1),1);
end

%Replacement of missing values
jtmd = repmat(1,1,size(x,2))-2*(sum(isnan(x))>0);
valmisdat = min(min(x))-0.5;
x(isnan(x)) = valmisdat;
valmd = repmat(valmisdat,1,size(x,2));

%Actual calculations
disv=daisyc(size(x,1),size(x,2),x,valmd,jtmd,vtype,metr);

%Putting things together
result = struct('disv',disv,'metric',metric,'number',size(x,1));


