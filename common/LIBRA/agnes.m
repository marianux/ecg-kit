function result = agnes(x,method,vtype,stdize,metric,banner)

%AGNES is a agglomerative clustering algorithm. It returns a hierarchy of clusters.
%
%The algorithm is fully described in:
%   Kaufman, L. and Rousseeuw, P.J. (1990),
%   "Finding groups in data: An introduction to cluster analysis",
%   Wiley-Interscience: New York (Series in Applied Probability and
%   Statistics), ISBN 0-471-87876-6.
%
% Required input arguments:
%        x : Data matrix (rows = observations, columns = variables)
%            or Dissimilarity matrix (if number of rows equals 1)
%   method : method used
%            Possible values are 1 : Group Average
%                                2 : Single Linkage
%                                3 : Complete Linkage
%                                4 : Ward
%                                5 : Weighted Average
%    vtype : Variable type vector (length equals number of variables)
%            Possible values are 1  Asymmetric binary variable (0/1)
%                                2  Nominal variable (includes symmetric binary)
%                                3  Ordinal variable
%                                4  Interval variable
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
%     banner : draws picture
%              Possible values are 0 : do not create a banner (default)
%                                  1 : create a banner
% I/O:
%   result=agnes(x,method,vtype,stdize,metric,banner)
%
% Example (subtracted from the referenced book)
%   load agricul.mat
%   result=agnes(agricul,1,[4 4],0,0,1);
%
% The output of AGNES is a structure containing:
%   result.x           : input matrix x
%   result.diss        : whether the inputmatrix x is a dissimilarity matrix
%                        or not
%   result.dys         : calculated dissimilarities (read row by row from the
%                        lower dissimilarity matrix, without the elements of
%                        the diagonal)
%   result.metric      : metric used 
%   result.stdize      : standardisation used
%   result.number      : number of observations
%   result.method      : method used
%   result.objectorder : order of objects
%   result.heights     : dissimilarity between the 2 relative clusters
%                        (=length of banner)
%   result.ac          : agglomerative coefficient
%   result.merge       : a (n-1) by 2 matrix related to the merge
%
% And AGNES will create the plot banner if banner equals 1.
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at:
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Wai Yan Kong 
% Created on 05/2006
% Last Revision: 19/09/2006


%Checking and filling out the inputs
if (nargin<2)
    error('Two input arguments required (data or dissimilarity vector and method)')
elseif ((nargin<3) & (size(x,1)~=1))
    error('Three input arguments required (data matrix x, method and vtype)')
    % so, only datamatrix x and method as input
elseif (nargin<3)
    metri ='unknown';
    metr='unknown';
    stdize = 0;
    banner = 0;
    % so, only dissim matrix x and method as input
elseif (nargin<4)
    stdize = 0;
    banner = 0;
    if (sum(vtype)~=4*size(x,2))
        metr ='mixed';
        metri='mixed';
    else
        metr ='eucli';
        metri='euclidean';
    end
    % so, only datamatrix or dissimilarity matrix x , method and vtype as
    % input
elseif (nargin<5)
    banner = 0;
    if (sum(vtype)~=4*size(x,2))
        metr ='mixed';
        metri='mixed';
    else
        metr ='eucli';
        metri='euclidean';
    end
    % so, only datamatrix or dissimilarity matrix x, method, vtype and
    % stdize as input
elseif (nargin<6)
    banner=0;
elseif (nargin>6)
    error('Too many input arguments')
end

% defining metric (for 4 input arguments) and diss
if (nargin>=5)
    if (metric==0)
        metr='eucli';
        metri='euclidean';
    elseif (metric==1)
        metr='manha';
        metri='manhattan';
    elseif (metric==2)
        metr='mixed';
        metri='mixed';
    else
        error('metric must be 0,1 or 2')
    end
end

if ((size(x,1)~=1))
    diss=0;
    dissi='x is no dissimilarity matrix';
else
    diss=1;
    dissi='x is a dissimilarity matrix';
end

%Standardization
if ((stdize==1) & (strcmp(metr,'eucli')| strcmp(metr,'manha')) & (diss==0))
    x = ((x - repmat(mean(x),size(x,1),1))./(repmat(std(x),size(x,1),1)));
    standardisation='standardisation by mean';
elseif ((stdize==2) & (strcmp(metr,'eucli')| strcmp(metr,'manha')) & (diss==0))
    x = ((x - repmat(median(x),size(x,1),1))./(repmat(mad(x),size(x,1),1)));
    standardisation='standardisation by median';
elseif(stdize==0)
    standardisation='no standardisation';
elseif (stdize==1 | stdize==2)
    standardisation='no standardisation(not enough num var or x is a diss matrix)';
elseif (nargin<=3)
    standardisation='no standardisation';
else
    error('stdize must be 0,1 or 2');
end

% defining dissimilarity matrix, number and method
if (diss==1)
    disv=x;
    number=(1+sqrt(1+8*size(x,2)))/2; %number of observations
     % checking for missing values in the dissimilarity matrix
    if any(isnan(disv))
        error('There are missing value(s) in the dissimilarity matrix!')
    end
    % checking the dimensions of the dissimilarity matrix
    if mod(number,fix(number))~=0
         error(['The dimension of the dissimilarity matrix is not correct!'])
    end
else
    resl=daisy(x,vtype,metr);
    disv=resl.disv;
    number=size(x,1);
end

if (method==1)
    namemethod='Group Average';
elseif (method==2)
    namemethod='Single Linkage';
elseif (method==3)
    namemethod='Complete Linkage';
elseif (method==4)
    namemethod='Ward';
elseif (method==5)
    namemethod='Weighted Average';
else
    error('method must be 1,2,3,4, or 5')
end

%Actual calculations
[ner,ban,coef,merge,dys]=twinsc(number,[0 disv]',method,1);

% We want ban to be a vector of length n-1
ban2=zeros(1,(number-1));
for i = 1:(number-1)
    ban2(i) = ban(i+1);
end

% We want merge to be a (n-1) by 2 matrix
merge2=ones(number-1,2);
for i = 1:(number-1)
    merge2(i,:) = merge(2*i-1:2*i);
end

% Create a banner
if (banner==1)
    Y=ban2;
    Y1=fliplr(Y);
    whitebg([0.4 0.5 0.75]);
    barh(Y1,1,'w');
    title('Banner of Agnes')
    xlabel('Height')
    YT=0.5:number;
    set(gca,'YTick',YT)
    set(gca,'YTickLabel',fliplr(ner))
    set(gca,'YAxisLocation','right')
    axis([min([ban2 0]),max([ban2 0]),0.5,number-0.5])
    set(gcf,'Color',[0.855388 0.855388 0.855388]);
elseif ((banner~=0) & (banner~=1) & (nargin==6))
    error('banner must equals 0 or 1')
end

%Putting things together
result = struct('x',x,'diss',dissi,'dys',dys,'metric',metri,...
    'stdize',standardisation,'number',number,...
    'method',namemethod,'objectorder',ner,'heights',ban2,...
    'ac',coef,'merge',merge2);
if diss
    result=rmfield(result, 'x');
end

