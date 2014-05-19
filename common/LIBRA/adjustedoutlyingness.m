function result = adjustedoutlyingness(x,varargin)

%ADJUSTEDOUTLYINGNESS calculates the 'Skewness-Adjusted Outlyingness'.
% The method searches for outliers in multivariate skewed data
% (thus without assuming elliptical symmetric).
% It is based on the outlyingness measure of Stahel and Donoho (1981, 1982)
% and on the skewness-adjusted boxplot of Hubert and Vandervieren (2007).
%
% The Skewness-Adjusted Outlyingness is described in:
%    Hubert, M., and Van der Veeken, S. (2008),
%    "Outlier detection for skewed data",
%    Journal of Chemometrics, 22, 235-246.
%
% The method can be useful as preprocessing in
% FASTICA (www.cis.hut.fi/projects/ica/fastica/), see
%
%    Brys, G., Hubert, M., and Rousseeuw, P.J. (2005),
%    "A Robustification of Independent Component Analysis",
%    Journal of Chemometrics, 19, 364-375.
%
%
% Required input arguments:
%    x : Data matrix (rows=observations, columns=variables)
%
% Optional input arguments:
%    ndir   : Number of directions in which de outlyingness will be computed
%             (default = 250*number of variables)
%   classic : If equal to one, the classical Stahel Donoho outlyingness is
%             calculated as well (default = 0)
%
% When the number of observations n is sufficiently large compared to the
% dimension p, i.e. when n > 5*p, directions are generated orthogonal to the 
% hyperplane through p random observations. The adjusted outlyingness is
% then affine invariant. When n is small compared to p, we take random
% directions through two data points. This procedure is orthogonal
% invariant.
%
% I/O:
%    result=adjustedoutlyingness(x,'ndir',250);
%
% Example:
%    x = [chi2rnd(2,1000,1) trnd(3,1000,1)];
%    x(1:10,:) = mvnrnd([-2 -3],eye(2)/10,10);
%    result = adjustedoutlyingness(x);
%
% The output of ADJUSTEDOUTLYINGNESS is a structure containing:
%    result.adjout      : skewness-adjusted outlyingness values for all observations
%    result.cutoff      : cutoff value for the AO-values
%    result.flag        : The observations whose AO-value exceeds the cutoff value can be
%                         considered as outliers and receive a flag equal to zero. The regular observations
%                         receive a flag 1.
%    result.classic     : If the input argument 'classic' is equal to one, this structure
%                         contains results of the classical analysis:
%                         .outl: Stahel-Donoho outlyingness for all observations
%                         .cutoff: cutoff value for the outlyingness values.
%                         .flag: The observations whose outlyingness exceeds the cutoff receive flag zero, the others receive flag 1.
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at:
%              http://wis.kuleuven.be/stat/robust
%
% Written by Guy Brys, Mia Hubert, Stephan Van der Veeken, Tim Verdonck
% Last Update: 09/06/2008

if (nargin<1)
    error('Input matrix x is required');
end
[n,p]=size(x);
counter=1;
default=struct('ndir',250*p,'a',-4,'b',3,'classic',0);
list=fieldnames(default);
result=default;
IN=length(list);
i=1;
%reading the user's input
if nargin>1
    %
    %placing inputfields in array of strings
    %
    for j=1:nargin-1
        if rem(j,2)~=0
            chklist{i}=varargin{j};
            i=i+1;
        end
    end
    %
    %Checking which default parameters have to be changed
    % and keep them in the structure 'result'.
    %
    while counter<=IN
        index=strmatch(list(counter,:),chklist,'exact');
        if ~isempty(index) %in case of similarity
            for j=1:nargin-1 %searching the index of the accompanying field
                if rem(j,2)~=0 %fieldnames are placed on odd index
                    if strcmp(chklist{index},varargin{j})
                        I=j;
                    end
                end
            end
            result=setfield(result,chklist{index},varargin{I+1});
            index=[];
        end
        counter=counter+1;
    end
end

ndir=result.ndir;
a=result.a;
b=result.b;
classic=result.classic;
if (ndir<=0)
    error('The number of directions should be positive.');
end

if isempty(ndir)
    ndir=250*p;
end

% a and b must be numeric scalars, classic should be a binary variable
if isempty(a)
    a = -4;
elseif ~isscalar(a) || ~isnumeric(a)
    error('The ''a'' parameter value must be a numeric scalar.');
end

if isempty(b)
    b = 3;
elseif ~isscalar(b) || ~isnumeric(b)
    error('The ''b'' parameter value must be a numeric scalar.');
end

if isempty(classic)
    classic = 0;
elseif ~isscalar(classic) || ~ismember(classic,0:1)
    error('Invalid value for ''classic'' parameter.');
end

if n>5*p
    B=[];
    for i =1:ndir
        xx = randperm(n);
        P = x(xx(1:p),:);
        if (rank(P)==p)
            B = [B ; (P\ones(p,1))'];
        end
    end
else
    B=twopoints(x,ndir,0);
end
for i=1:size(B,1)
    Bnorm(i)=norm(B(i,:),2);
end
Bnormr=Bnorm(Bnorm > 1.e-12);
B=B(Bnorm > 1.e-12,:);
A=diag(1./Bnormr)*B;

%Looking in ndir directions for skewness-adjusted outlyingness
Y=x*A';
[n,p]=size(Y);
YZ=zeros(n,p);

tmc = mc(Y);
h=find(abs(tmc)==1);
if(sum(h)>1)
    error('There are too many ties in the data. The adjusted outlyingness can not be computed.')
end
tme = median(Y);
tiq = iqr(Y);

s=find(tiq==0);
if(sum(s)>1)
    error('There are too many ties in the data. The adjusted outlyingness can not be computed.')
end

tp1 = prctile(Y,25);
tp3 = prctile(Y,75);
if (tmc>=0)
    tup = (tp3+1.5*exp(b*tmc).*tiq)-tme;%3
    tlo = tme-(tp1-1.5*exp(a*tmc).*tiq);%-4
else
    tup = (tp3+1.5*exp(-a*tmc).*tiq)-tme;%4
    tlo = tme-(tp1-1.5*exp(-b*tmc).*tiq);%-3
end
for k = 1:p
    ttup = sort(-Y(Y(:,k)<(tup(k)+tme(k)),k));
    tup(k) = -ttup(1)-tme(k);
    ttlo = sort(Y(Y(:,k)>(tme(k)-tlo(k)),k));
    tlo(k) = tme(k)-ttlo(1);
end
tmp = ((Y>=repmat(tme,n,1)).*(repmat(tup,n,1)) + (Y<repmat(tme,n,1)).*(repmat(tlo,n,1)));
YZ = abs((Y-repmat(tme,n,1)))./tmp;
adjout = max(YZ,[],2);

mcadjout=mc(adjout);
if mcadjout>0
    cutoff = prctile(adjout,75)+1.5*exp(b*mcadjout)*iqr(adjout);
else
    cutoff = prctile(adjout,75)+1.5*iqr(adjout);
end
ttup=sort(-adjout(adjout<cutoff));
cutoff=-ttup(1);
flag=(adjout<=cutoff);

if classic==1
    % The classical estimates are computed.
    umad=madc(Y);
    g=find(umad==0);
    if(sum(g)>1)
        error('There are too many ties in the data. The Stahel-Donoho outlyingness can not be computed.')
    end

    clYZ=abs((Y-repmat(tme,n,1)))./repmat(umad,n,1);
    classout=max(clYZ,[],2);

    %classical chi-square cutoff and flag
    quant=chi2inv(0.99,size(x,2));
    classcutoff=sqrt(quant);
    classflag=(classout<=classcutoff);

    classic=struct('outl',classout,'cutoff',classcutoff,'flag',classflag);

end

%Putting things together
result = struct('adjout',adjout,'cutoff',cutoff,'flag',flag,'classic',classic);