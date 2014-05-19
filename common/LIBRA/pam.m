function result = pam(x,kclus,vtype,stdize,metric,silhplot)

%PAM is the Partitioning Around Medoids clustering algorithm.
% It returns a list representing a clustering of the data into kclus
% clusters based on the search for kclus representative objects or medoids among the observations of
% the data set.
%
% The algorithm is fully described in:
%   Kaufman, L. and Rousseeuw, P.J. (1990),
%   "Finding groups in data: An introduction to cluster analysis",
%   Wiley-Interscience: New York (Series in Applied Probability and
%   Statistics), ISBN 0-471-87876-6.
%
% Required input arguments:
%       x : Data matrix (rows = observations, columns = variables)
%           or Dissimilarity matrix (if number of columns equals 1)
%   kclus : The number of desired clusters
%   vtype : Variable type vector (length equals number of variables)
%           Possible values are 1  Asymmetric binary variable (0/1)
%                               2  Nominal variable (includes symmetric binary)
%                               3  Ordinal variable
%                               4  Interval variable
%          (if x is a dissimilarity matrix vtype is not required.)
%
% Optional input arguments:
%     stdize : standardise the variables given by the x-matrix
%              Possible values are 0 : no standardisation (default)
%                                  1 : standardisation by the mean
%                                  2 : standardisation by the median
%              (if x is a dissimilarity matrix, stdize is ignored)
%     metric : Metric to be used 
%              Possible values are 'eucli' Euclidian (all interval variables, default)
%                                  'manha' Manhattan
%                                  'mixed' Mixed (not all interval variables, default)
%              (if x is a dissimilarity matrix, metric is ignored)
%   silhplot : draws picture
%              Possible values are 0 : do not create a silhouette plot (default)
%                                  1 : create a silhouette plot
%
% I/O:
%   result=pam(x,kclus,vtype,'eucli',silhplot)
%
% Example (subtracted from the referenced book)
%   load ruspini.mat
%   result=pam(ruspini,2,[4 4],1);
%
% The output of PAM is a structure containing:
%   result.dys        : dissimilarities (read row by row from the
%                       lower dissimilarity matrix)
%   result.metric     : metric used
%   result.number     : number of observations
%   result.ttd        : Average silhouette width per cluster
%   result.ttsyl      : Average silhouette width for dataset
%   result.idmed      : Id of medoid observations
%   result.obj        : Objective function at the first two iterations
%   result.ncluv      : Cluster membership for each observation
%   result.clusinf    : Matrix, each row gives numerical information for
%                       one cluster. These are the cardinality of the cluster
%                       (number of observations), the maximal and average
%                       dissimilarity between the observations in the cluster
%                       and the cluster's medoid, the diameter of the cluster
%                       (maximal dissimilarity between two observations of the
%                       cluster), and the separation of the cluster (minimal
%                       dissimilarity between an observation of the cluster
%                       and an observation of another cluster).
%   result.sylinf     : Matrix, with for each observation i the cluster to
%                       which i belongs, as well as the neighbor cluster of i
%                       (the cluster, not containing i, for which the average
%                       dissimilarity between its observations and i is minimal),
%                       and the silhouette width of i.
%   result.nisol      : Vector, with for each cluster specifying whether it is
%                       an isolated cluster (L- or L*-clusters) or not isolated.
%                       A cluster is an L*-cluster iff its diameter is smaller than
%                       its separation.  A cluster is an L-cluster iff for each
%                       observation i the maximal dissimilarity between i and any
%                       other observation of the cluster is smaller than the minimal
%                       dissimilarity between i and any observation of another cluster.
%                       Clearly each L*-cluster is also an L-cluster.
%
% And PAM will create the silhouette plot if silhplot equals 1 (an empty bar indicated by
%                       zero is a sparse between two clusters).
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at:
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Guy Brys (May 2006)

%Checking and filling in the inputs
res1=[];
if (nargin<2)
    error('Two input arguments required')
elseif ((nargin<3) & (size(x,2)~=1) & (size(x,1)~=1))
    error('Three input arguments required')
elseif (nargin<3)
    if (size(x,2)==1)
        x = x';
    end
    res1.metric = 'unknown';
    res1.disv = x;
    lookup=seekN(x);
    res1.number = lookup.numb;   %(1+sqrt(1+8*size(x,1)))/2;
    stdize = 0;
    silhplot = 0;
elseif (nargin<4)
    stdize = 0;
    silhplot = 0;
    if (sum(vtype)~=4*size(x,2))
        metric = 'mixed';
    else
        metric = 'eucli';
    end
elseif (nargin<5)
    silhplot = 0;
    if (sum(vtype)~=4*size(x,2))
        metric = 'mixed';
    else
        metric = 'eucli';
    end
elseif (nargin<6)
    silhplot = 0;
end

%Replacement of missing values
for i=1:size(x,1)
    A=find(isnan(x(i,:)));
    if (~(isempty(A)))
        for j=A
            valmisdat=0;
            for c=1:size(x,2)
                if (c~=j)
                    [a,b] = sort(x(:,c));
                    if (~isempty(b(find(a==x(i,c)))))
                        valmisdat=valmisdat+find(a==x(i,c));
                    end
                end
            end
            x(i,j)=prctile(x(isnan(x(:,j))==0,j),100*valmisdat/(size(x,1)*(size(x,2)-1)));
        end
    end
end

%Standardization
if ((stdize==1) & (strcmp(metric,'eucli')))
    x = ((x - repmat(mean(x),size(x,1),1))./(repmat(std(x),size(x,1),1)));
elseif ((stdize==2) & (strcmp(metric,'eucli')))
    x = ((x - repmat(median(x),size(x,1),1))./(repmat(mad(x),size(x,1),1)));
end

%Calculating the dissimilarities with daisy
if (isempty(res1))
    res1=daisy(x,vtype,metric);
end

%Actual calculations (the second for latter use with CLUSPLOT)
[dys,ttd,ttsyl,idmed,obj,ncluv,clusinf,sylinf,nisol]=pamc(res1.number,kclus,[0 res1.disv]');
dys=res1.disv(lowertouppertrinds(res1.number));

%Create a silhouetteplot
if (silhplot==1)
    Y=sylinf(:,3);
    Y1=flipdim(Y,1);
    whitebg([1 1 1]);
    % we calculate b="a but with a bar with length zero if the objects
    % are from another cluster"
    % and h="objects but with a 0 between 2 clusters"="g with a 0 if
    % it is a sparse between 2 clusters"
    a=flipdim(Y1,1);
    b=[];
    g=sylinf(:,4);
    f=sylinf(:,1)-1;
    for j=1:res1.number
        b(j+f(j))=a(j);
        h(j+f(j))=g(j);
    end
    b1=flipdim(b,2);
    h1=flipdim(h,2);
    % we use this b1 and h1 to plot the barh (instead of a and g)
    barh(b1,1);
    title 'Silhouette Plot of Pam' ;
    xlabel('Silhouette width');
    YT=1:res1.number+(sylinf(res1.number,1)-1);
    set(gca,'YTick',YT);
    set(gca,'YTickLabel',h1);
    axis([min([Y' 0]),max([Y' 0]),0.5,res1.number+0.5+f(res1.number)]);
elseif ((silhplot~=0) & (silhplot~=1) & (nargin==6))
    error('silhplot must equals 0 or 1')
end


%Putting things together
result = struct('dys',dys,'metric',res1.metric,'number',res1.number,...
    'ttd',ttd,'ttsyl',ttsyl,'idmed',idmed,'obj',obj,'ncluv',ncluv,...
    'clusinf',clusinf,'sylinf',sylinf,'nisol',nisol,'x',x);

%------------
%SUBFUNCTIONS

function dv = lowertouppertrinds(n)

dv=[];
for i=0:(n-2)
    dv = [dv cumsum(i:(n-2))+repmat(1+sum(0:i),1,n-i-1)];
end

%---
function outn = seekN(x)

ok=0;
numb=0;
k=size(x,2);
sums=cumsum(1:k);
for i=1:k
    if(sums(i)==k)

        numb=i+1;
        ok=1;
    end
end
outn=struct('numb',numb,'ok',ok);
