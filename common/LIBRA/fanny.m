function result = fanny(x,kclus,vtype,metric,silhplot)

%FANNY is a fuzzy clustering algorithm. It returns a list representing a fuzzy clustering of the data
% into kclus clusters.
%
%The algorithm is fully described in:
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
%   result=fanny(x,kclus,vtype,'eucli',silhplot)
%
% Example (subtracted from the referenced book)
%   load country.mat
%   result=fanny(country,2,[4 4]);
%
% The output of FANNY is a structure containing:
%   result.dys        : dissimilarities (read row by row from the
%                       lower dissimilarity matrix)
%   result.metric     : metric used 
%   result.number     : number of observations
%   result.pp         : Membership coefficients for each observation
%   result.coeff      : Dunn's partition coefficient (and normalized version)
%   result.ncluv      : A vector with length equal to the number of observations,
%                       giving for each observation the number of the cluster to
%                       which it has the largest membership
%   result.obj        : Objective function and the number of iterations the
%                       fanny algorithm needed to reach this minimal value
%   result.sylinf     : Matrix, with for each observation i the cluster to
%                       which i belongs, as well as the neighbor cluster of i
%                       (the cluster, not containing i, for which the average
%                       dissimilarity between its observations and i is minimal),
%                       and the silhouette width of i.
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at:
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Guy Brys and Wai Yan Kong (May 2006)

%Checking and filling in the inputs
res1=[];
if (nargin<2)
    error('Two input arguments required')
elseif ((nargin<3) & (size(x,2)~=1))
    error('Three input arguments required')
elseif (nargin<3)
    res1.metric = 'unknown';
    res1.disv = x';
    lookup=seekN(x);
    res1.number = lookup.numb; %(1+sqrt(1+8*size(x,1)))/2;
    silhplot = 0;
elseif (nargin<4)
    silhplot = 0;
    if (sum(vtype)~=4*size(x,2))
        metric = 'mixed';
    else
        metric = 'eucli';
    end
elseif (nargin<5)
    silhplot = 0;
end

%Calculating the dissimilarities with daisy
%For fanny the second command is also required
if (isempty(res1))
    res1=daisy(x,vtype,metric);
end
res1.disv=res1.disv(lowertouppertrinds(res1.number));

%Actual calculations
[pp,coeff,clu,obj,sylinf]=fannyc(res1.number,kclus,[0 res1.disv]');
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
    title 'Silhouette Plot of Fanny' ;
    xlabel('Silhouette width');
    YT=1:res1.number+(sylinf(res1.number,1)-1);
    set(gca,'YTick',YT);
    set(gca,'YTickLabel',h1);
    axis([min([Y' 0]),max([Y' 0]),0.5,res1.number+0.5+f(res1.number)]);
elseif ((silhplot~=0) & (silhplot~=1) & (nargin==5))
    error('silhplot must equals 0 or 1')
end

%Putting things together
result = struct('dys',res1.disv,'metric',res1.metric,...
    'number',res1.number,'pp',pp,...
    'coeff',coeff,'ncluv',clu,'obj',obj,'sylinf',sylinf);

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
k=size(x,1);
sums=cumsum(1:k);
for i=1:k
    if(sums(i)==k)
        numb=i+1;
        ok=1;
    end
end
outn=struct('numb',numb,'ok',ok);
