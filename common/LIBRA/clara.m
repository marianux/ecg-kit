function result = clara(x,kclus,vtype,stdize,metric,nsamp,sampsize)

%CLARA is the 'Clustering Large Applications' clustering algorithm.
% It returns a list representing a clustering of the data
% into kclus clusters following the clara algorithm which is
% designed for large data sets.
%
%The algorithm is fully described in:
%   Kaufman, L. and Rousseeuw, P.J. (1990),
%   "Finding groups in data: An introduction to cluster analysis",
%   Wiley-Interscience: New York (Series in Applied Probability and
%   Statistics), ISBN 0-471-87876-6.
%
% Required input arguments:
%       x : Data matrix (rows = observations, columns = variables)
%   kclus : The number of desired clusters
%   vtype : Variable type vector (length equals number of variables)
%           Possible values are 1  Asymmetric binary variable (0/1)
%                               2  Nominal variable (includes symmetric binary)
%                               3  Ordinal variable
%                               4  Interval variable
%
% Optional input arguments:
%     stdize : standardise the variables given by the x-matrix
%              Possible values are 0 : no standardisation (default)
%                                  1 : standardisation by the mean
%                                  2 : standardisation by the median
%     metric : Metric to be used 
%              Possible values are 'eucli' Euclidian (all interval variables, default)
%                                  'manha' Manhattan
%                                  'mixed' Mixed (not all interval variables, default)
%      nsamp : Number of samples to be drawn from the data set
%   sampsize : Number of observations in each sample (should be higher
%              than the number of clusters and lower than the number of
%              observations)
%
% I/O:
%   result=clara(x,kclus,vtype,'eucli',5,40+2*kclus)
%
% Example (subtracted from the referenced book)
%   load obj200.mat
%   result=clara(obj200,3,[4 4]);
%
% The output of CLARA is a structure containing:
%   result.dysobs     : dissimilarities for each observation with the medoids
%   result.metric     : metric used
%   result.number     : number of observations
%   result.idmed      : Id of medoid observations
%   result.ncluv      : A vector with length equal to the number of observations,
%                       giving for each observation the number of the cluster to
%                       which it belongs
%   result.obj        : Objective function for the best subsample
%   result.clusinf    : Matrix, each row gives numerical information for
%                       one cluster. These are the cardinality of the cluster
%                       (number of observations), the maximal and average
%                       dissimilarity between the observations in the cluster
%                       and the cluster's medoid, the diameter of the cluster
%                       (maximal dissimilarity between two observations of the
%                       cluster), and the separation of the cluster (minimal
%                       dissimilarity between an observation of the cluster
%                       and an observation of another cluster).
%   result.sylinf     : Matrix based on the best subsample, with for each
%                       observation i of this subsample the cluster to
%                       which i belongs, as well as the neighbor cluster of i
%                       (the cluster, not containing i, for which the average
%                       dissimilarity between its observations and i is minimal),
%                       and the silhouette width of i.
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at:
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Guy Brys (May 2006)

%Checking and filling out the inputs
if (nargin<3)
    error('Three input arguments required')
elseif (nargin<4)
    stdize = 0;
    if (sum(vtype)~=4*size(x,2))
        metric = 'mixed';
    else
        metric = 'eucli';
    end
    nsamp=5;
    sampsize=40+2*kclus;
elseif (nargin<5)
    if (sum(vtype)~=4*size(x,2))
        metric = 'mixed';
    else
        metric = 'eucli';
    end
    nsamp=5;
    sampsize=40+2*kclus;
elseif (nargin<6)
    nsamp=5;
    sampsize=40+2*kclus;
elseif (nargin<7)
    sampsize=40+2*kclus;
end

%Standardization
if (stdize==1)
    x = ((x - repmat(mean(x),size(x,1),1))./(repmat(std(x),size(x,1),1)));
elseif (stdize==2)
    x = ((x - repmat(median(x),size(x,1),1))./(repmat(mad(x),size(x,1),1)));
end

%Actual calculations
obj = Inf;
for i=1:nsamp
    sampindex = randperm(size(x,1));
    restemp = pam(x(sampindex(1:sampsize),:),kclus,vtype,metric);
    if (restemp.obj(1)<obj)
        obj = restemp.obj(1);
        idmed = sampindex(restemp.idmed);
    end
end

%Calculating some extra dissimilarities for output
for i=1:size(x,1)
    for j=1:kclus
        distemp = daisy(x([i idmed(j)],:),vtype,metric);
        disv(i,j) = distemp.disv(1);
    end
    [zz,clu(i)] = min(disv(i,:));
end
mindisv=[];
for j=1:kclus
    clusinf(j,1) = length(find(clu==j));
    clusinf(j,2) = max(disv(find(clu==j),j));
    clusinf(j,3) = mean(disv(find(clu==j),j));
end
for i=1:kclus
    for j=1:kclus
        distemp = daisy(x([idmed(i) idmed(j)],:),vtype,metric);
        mindisv = [mindisv distemp.disv(1)];
    end
end
clusinf(:,4) = clusinf(:,2)/min(mindisv(mindisv~=0));

%Putting things together
result = struct('dysobs',disv,'metric',metric,'number',size(x,1),...
    'idmed',idmed,'ncluv',clu,'obj',obj,'clusinf',clusinf,...
    'sylinf',restemp.sylinf,'x',x);







