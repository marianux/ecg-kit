function [rew,raw]=mcdcov(x,varargin)

%MCDCOV computes the MCD estimator of a multivariate data set.  This
% estimator is given by the subset of h observations with smallest covariance
% determinant.  The MCD location estimate is then the mean of those h points,
% and the MCD scatter estimate is their covariance matrix.  The default value
% of h is roughly 0.75n (where n is the total number of observations), but the
% user may choose each value between n/2 and n. Based on the raw estimates,
% weights are assigned to the observations such that outliers get zero weight.
% The reweighted MCD estimator is then given by the mean and covariance matrix
% of the cases with non-zero weight. To compute the MCD estimator,
% the FASTMCD algorithm is used.
%
% The MCD method is intended for continuous variables, and assumes that
% the number of observations n is at least 5 times the number of variables p.
% If p is too large relative to n, it would be better to first reduce
% p by variable selection or robust principal components (see the functions
% robpca.m and rapca.m).
%
% The MCD method was introduced in:
%
%   Rousseeuw, P.J. (1984), "Least Median of Squares Regression,"
%   Journal of the American Statistical Association, Vol. 79, pp. 871-881.
%
% The MCD is a robust method in the sense that the estimates are not unduly
% influenced by outliers in the data, even if there are many outliers.
% Due to the MCD's robustness, we can detect outliers by their large
% robust distances. The latter are defined like the usual Mahalanobis
% distance, but based on the MCD location estimate and scatter matrix
% (instead of the nonrobust sample mean and covariance matrix).
%
% The FASTMCD algorithm uses several time-saving techniques which
% make it available as a routine tool to analyze data sets with large n,
% and to detect deviating substructures in them. A full description of the
% algorithm can be found in:
%
%   Rousseeuw, P.J. and Van Driessen, K. (1999), "A Fast Algorithm for the
%   Minimum Covariance Determinant Estimator," Technometrics, 41, pp. 212-223.
%
% An important feature of the FASTMCD algorithm is that it allows for exact
% fit situations, i.e. when more than h observations lie on a (hyper)plane.
% Then the program still yields the MCD location and scatter matrix, the latter
% being singular (as it should be), as well as the equation of the hyperplane.
%
%
% Required input argument:
%    x : a vector or matrix whose columns represent variables, and rows represent observations.
%        Missing values (NaN's) and infinite values (Inf's) are allowed, since observations (rows)
%        with missing or infinite values will automatically be excluded from the computations.
%
% Optional input arguments:
%       cor : If non-zero, the robust correlation matrix will be
%             returned. The default value is 0.
%         h : The quantile of observations whose covariance determinant will
%             be minimized.  Any value between n/2 and n may be specified.
%             The default value is 0.75*n.
%     alpha : (1-alpha) measures the fraction of outliers the algorithm should
%             resist. Any value between 0.5 and 1 may be specified. (default = 0.75)
%    ntrial : The number of random trial subsamples that are drawn for
%             large datasets. The default is 500.
%     plots : If equal to one, a menu is shown which allows to draw several plots,
%             such as a distance-distance plot. (default)
%             If 'plots' is equal to zero, all plots are suppressed.
%             See also makeplot.m
%   classic : If equal to one, the classical mean and covariance matrix are computed as well.
%             (default = 0)
%
% Input arguments for advanced users:
%     Hsets : Instead of random trial h-subsets (default, Hsets = []), Hsets makes it possible to give certain
%             h-subsets as input. Hsets is a matrix that contains the indices of the observations of one
%             h-subset as a row.
%    factor : If not equal to 0 (default), the consistency factor is adapted. Only usefull in case of the
%             kmax approach.
%
% I/O: result=mcdcov(x,'alpha',0.75,'h',h,'ntrial',500)
%  If only one output argument is listed, only the final result ('result')
%  is returned.
%  The user should only give the input arguments that have to change their default value.
%  The name of the input arguments needs to be followed by their value.
%  The order of the input arguments is of no importance.
%
% Examples: [rew,raw]=mcdcov(x);
%           result=mcdcov(x,'h',20,'plots',0);
%           [rew,raw]=mcdcov(x,'alpha',0.8,'cor',0)
%
% The output structure 'raw' contains intermediate results, with the following
% fields :
%
%     raw.center : The raw MCD location of the data.
%        raw.cov : The raw MCD covariance matrix (multiplied by a consistency factor).
%        raw.cor : The raw MCD correlation matrix, if input argument 'cor' was non-zero.
%         raw.wt : Weights based on the estimated raw covariance matrix 'raw.cov' and
%                  the estimated raw location 'raw.center' of the data. These weights determine
%                  which observations are used to compute the final MCD estimates.
%  raw.objective : The determinant of the raw MCD covariance matrix.
%
% The output structure 'rew' contains the final results, namely:
%
%       rew.center : The robust location of the data, obtained after reweighting, if
%                    the raw MCD is not singular.  Otherwise the raw MCD center is
%                    given here.
%          rew.cov : The robust covariance matrix, obtained after reweighting, if the raw MCD
%                    is not singular.  Otherwise the raw MCD covariance matrix is given here.
%          rew.cor : The robust correlation matrix, obtained after reweighting, if
%                    options.cor was non-zero.
%            rew.h : The number of observations that have determined the MCD estimator,
%                    i.e. the value of h.
%    rew. Hsubsets : A structure that contains Hopt and Hfreq:
%                    Hopt  : The subset of h points whose covariance matrix has minimal determinant,
%                            ordered following increasing robust distances.
%                    Hfreq : The subset of h points which are the most frequently selected during the whole
%                            algorithm.
%        rew.alpha : (1-alpha) measures the fraction of outliers the algorithm should
%                    resist.
%           rew.md : The distance of each observation from the classical
%                    center of the data, relative to the classical shape
%                    of the data. Often, outlying points fail to have a
%                    large Mahalanobis distance because of the masking
%                    effect.
%           rew.rd : The distance of each observation to the final,
%                    reweighted MCD center of the data, relative to the
%                    reweighted MCD scatter of the data.  These distances allow
%                    us to easily identify the outliers. If the reweighted MCD
%                    is singular, raw.rd is given here.
%       rew.cutoff : Cutoff values for the robust and mahalanobis distances
%         rew.flag : Flags based on the reweighted covariance matrix and the
%                    reweighted location of the data.  These flags determine which
%                    observations can be considered as outliers. If the reweighted
%                    MCD is singular, raw.wt is given here.
%       rew.method : A character string containing information about the method and
%                    about singular subsamples (if any).
%        rew.plane : In case of an exact fit, rew.plane contains the coefficients
%                    of a (hyper)plane a_1(x_i1-m_1)+...+a_p(x_ip-m_p)=0
%                    containing at least h observations, where (m_1,...,m_p)
%                    is the MCD location of these observations.
%      rew.classic : If the input argument 'classic' is equal to one, this structure
%                    contains results of the classical analysis: center (sample mean),
%                    cov (sample covariance matrix), md (Mahalanobis distances), class ('COV').
%        rew.class : 'MCDCOV'
%            rew.X : If x is bivariate, same as the x in the call to mcdcov,
%                    without rows containing missing or infinite values.
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at:
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Katrien Van Driessen and Bjorn Rombouts
% Revisions by Sanne Engelen, Sabine Verboven
% Last Update: 09/04/2004, 01/08/2007

% The FASTMCD algorithm works as follows:
%
%       The dataset contains n cases and p variables.
%       When n < 2*nmini (see below), the algorithm analyzes the dataset as a whole.
%       When n >= 2*nmini (see below), the algorithm uses several subdatasets.
%
%       When the dataset is analyzed as a whole, a trial subsample of p+1 cases
%       is taken, of which the mean and covariance matrix are calculated.
%       The h cases with smallest relative distances are used to calculate
%       the next mean and covariance matrix, and this cycle is repeated csteps1
%       times. For small n we consider all subsets of p+1 out of n, otherwise
%       the algorithm draws 500 random subsets by default.
%       Afterwards, the 10 best solutions (means and corresponding covariance
%       matrices) are used as starting values for the final iterations.
%       These iterations stop when two subsequent determinants become equal.
%       (At most csteps3 iteration steps are taken.) The solution with smallest
%       determinant is retained.
%
%       When the dataset contains more than 2*nmini cases, the algorithm does part
%       of the calculations on (at most) maxgroup nonoverlapping subdatasets, of
%       (roughly) maxobs cases.
%
%       Stage 1: For each trial subsample in each subdataset, csteps1 (see below) iterations are
%       carried out in that subdataset. For each subdataset, the 10 best solutions are
%       stored.
%
%       Stage 2 considers the union of the subdatasets, called the merged set.
%       (If n is large, the merged set is a proper subset of the entire dataset.)
%       In this merged set, each of the 'best solutions' of stage 1 are used as starting
%       values for csteps2 (sse below) iterations. Also here, the 10 best solutions are stored.
%
%       Stage 3 depends on n, the total number of cases in the dataset.
%       If n <= 5000, all 10 preliminary solutions are iterated.
%       If n > 5000, only the best preliminary solution is iterated.
%       The number of iterations decreases to 1 according to n*p (If n*p <= 100,000 we
%       iterate csteps3 (sse below) times, whereas for n*p > 1,000,000 we take only one iteration step).
%

if rem(nargin-1,2)~=0
    error('The number of input arguments should be odd!');
end
% Assigning some input parameters
data = x;
raw.cor = [];
rew.cor = [];
rew.plane = [];
% The maximum value for n (= number of observations) is:
% nmax=50000;
nmax=realmax;
% The maximum value for p (= number of variables) is:
pmax=50;
% To change the number of subdatasets and their size, the values of
% maxgroup and nmini can be changed.
maxgroup=5;
nmini=300;
% The number of iteration steps in stages 1,2 and 3 can be changed
% by adapting the parameters csteps1, csteps2, and csteps3.
csteps1=2;
csteps2=2;
csteps3=100;
% dtrial : number of subsamples if not all (p+1)-subsets will be considered.
dtrial=500;

if size(data,1)==1
    data=data';
end

% Observations with missing or infinite values are ommitted.
ok=all(isfinite(data),2);
data=data(ok,:);
xx=data;
[n,p]=size(data);
% Some checks are now performed.
if n==0
    error('All observations have missing or infinite values.')
end
if n > nmax
    error(['The program allows for at most ' int2str(nmax) ' observations.'])
end
if p > pmax
    error(['The program allows for at most ' int2str(pmax) ' variables.'])
end
if n < p
    error('Need at least (number of variables) observations.')
end

%internal variables
hmin=quanf(0.5,n,p);
%Assiging default values
h=quanf(0.75,n,p);
default=struct('alpha',0.75,'h',h,'plots',1,'ntrial',dtrial,'cor',0,'seed',0,'classic',0,'Hsets',[],'factor',0);
list=fieldnames(default);
options=default;
IN=length(list);
i=1;
counter=1;

%Reading optional inputarguments
if nargin>2
    %
    % placing inputfields in array of strings
    %
    for j=1:nargin-1
        if rem(j,2)~=0
            chklist{i}=varargin{j};
            i=i+1;
        end
    end
    dummy=sum(strcmp(chklist,'h')+2*strcmp(chklist,'alpha'));
    switch dummy
        case 0 % defaultvalues should be taken
            alfa=options.alpha;
            h=options.h;
        case 3
            error('Both input arguments alpha and h are provided. Only one is required.')
    end
    %
    % Checking which default parameters have to be changed
    % and keep them in the structure 'options'.
    %
    while counter<=IN
        index=strmatch(list(counter,:),chklist,'exact');
        if ~isempty(index) % in case of similarity
            for j=1:nargin-2 % searching the index of the accompanying field
                if rem(j,2)~=0 % fieldnames are placed on odd index
                    if strcmp(chklist{index},varargin{j})
                        I=j;
                    end
                end
            end
            options=setfield(options,chklist{index},varargin{I+1});
            index=[];
        end
        counter=counter+1;
    end
    if dummy==1% checking inputvariable h
        % hmin is the minimum number of observations whose covariance determinant
        % will be minimized.
        if isempty(options.Hsets)
            if options.h < hmin
                disp(['Warning: The MCD must cover at least ' int2str(hmin) ' observations.'])
                disp(['The value of h is set equal to ' int2str(hmin)])
                options.h = hmin;
            elseif options.h > n
                error('h is greater than the number of non-missings and non-infinites.')
            elseif options.h < p
                error(['h should be larger than the dimension ' int2str(p) '.'])
            end
        end
        options.alpha=options.h/n;
    elseif dummy==2
        if options.alpha < 0.5
            options.alpha=0.5;
            mess=sprintf(['Attention (mcdcov.m): Alpha should be larger than 0.5. \n',...
                'It is set to 0.5.']);
            disp(mess)
        end
        if options.alpha > 1
            options.alpha=0.75;
            mess=sprintf(['Attention (mcdcov.m): Alpha should be smaller than 1.\n',...
                'It is set to 0.75.']);
            disp(mess)
        end
        options.h=quanf(options.alpha,n,p);
    end
end

h=options.h;  %number of regular datapoints on which estimates are based. h=[alpha*n]
plots=options.plots; %relevant plots are plotted
alfa=options.alpha; %percentage of regular observations
ntrial=options.ntrial; %number of subsets to be taken in the first step
cor=options.cor; %correlation matrix
seed=options.seed; %seed of the random generator
cutoff.rd=sqrt(chi2inv(0.975,p)); %cutoff value for the robust distance
cutoff.md=cutoff.rd; %cutoff value for the mahalanobis distance
Hsets = options.Hsets;
if ~isempty(Hsets)
    Hsets_ind = 1;
else
    Hsets_ind = 0;
end
factor = options.factor;
if factor == 0
    factor_ind = 0;
else
    factor_ind = 1;
end

% Some initializations.
rew.flag=repmat(NaN,1,length(ok));
raw.wt=repmat(NaN,1,length(ok));
raw.rd=repmat(NaN,1,length(ok));
rew.rd=repmat(NaN,1,length(ok));
rew.mahalanobis=repmat(NaN,1,length(ok));
rew.method=sprintf('\nMinimum Covariance Determinant Estimator.');
correl=NaN;

%    z    : if at least h observations lie on a hyperplane, then z contains the
%           coefficients of that plane.
% weights : weights of the observations that are not excluded from the computations.
%           These are the observations that don't contain missing or infinite values.
% bestobj : best objective value found.
z(1:p)=0;
weights=zeros(1,n);
bestobj=inf;

% The classical estimates are computed.
[Pcl,Tcl,Lcl,rcl,centerXcl,cXcl] = classSVD(data);
clasmean=cXcl;
clascov=Pcl*diag(Lcl)*Pcl';

if p < 5
    eps=1e-12;
elseif p <= 8
    eps=1e-14;
else
    eps=1e-16;
end

% The standardization of the data will now be performed.
med=median(data);
mad=sort(abs(data-repmat(med,n,1)));
mad=mad(h,:);
ii=min(find(mad < eps));
if length(ii)
    % The h-th order statistic is zero for the ii-th variable. The array plane contains
    % all the observations which have the same value for the ii-th variable.
    plane=find(abs(data(:,ii)-med(ii)) < eps)';
    meanplane=mean(data(plane,:));
    weights(plane)=1;
    if p==1
        rew.flag=weights;
        raw.wt=weights;
        [raw.center,rew.center]=deal(meanplane);
        [raw.cov,rew.cov,raw.objective]=deal(0);
        if plots
            rew.method=sprintf('\nUnivariate location and scale estimation.');
            rew.method=strvcat(rew.method,sprintf('%g of the %g observations are identical.',length(plane),n));
        end
    else
        z(ii)=1;
        rew.plane=z;
        covplane=cov(data(plane,:));
        [raw.center,raw.cov,rew.center,rew.cov,raw.objective,raw.wt,rew.flag,...
            rew.method]=displ(3,length(plane),weights,n,p,meanplane,covplane,rew.method,z,ok,...
            raw.wt,rew.flag,0,NaN,h,ii);
    end
    rew.Hsubsets.Hopt = plane;
    rew.Hsubsets.Hfreq = plane;
    %classical analysis?
    if options.classic==1
        classic.cov=clascov;
        classic.center=clasmean;
        classic.class='COV';
    else
        classic=0;
    end
    %assigning the output
    rewo=rew;rawo=raw;
    rew=struct('center',{rewo.center},'cov',{rewo.cov},'cor',{cor},'h',{h},'Hsubsets',{rewo.Hsubsets},'alpha',{alfa},...
        'flag',{rewo.flag},'plane', {rewo.plane},'method',{rewo.method},'class',{'MCDCOV'},'classic',{classic},'X',{xx});
    raw=struct('center',{rawo.center},'cov',{rawo.cov},'cor',{rawo.cor},'objective',{rawo.objective},...
        'wt',{rawo.wt},'class',{'MCDCOV'},'classic',{classic},'X',{x});
    if size(data,2)~=2
        rew=rmfield(rew,'X');
        raw=rmfield(raw,'X');
    end
    return
end
data=(data-repmat(med,n,1))./repmat(mad,n,1);

% The standardized classical estimates are now computed.
clmean=mean(data);
clcov=cov(data);

% The univariate non-classical case is now handled.
if p==1 & h~=n
    [rew.center, rewsca, weights,raw.center,raw.cov,rawdist,raw.wt,Hopt]=unimcd(data,h);
    rew.Hsubsets.Hopt = Hopt';
    rew.Hsubsets.Hfreq = Hopt';
    raw.rd=sqrt(rawdist');
    rew.cov=rewsca^2;
    raw.objective=raw.cov*prod(mad)^2;
    mah=(data-rew.center).^2/rew.cov;
    rew.rd=sqrt(mah');
    rew.flag= rew.rd <= cutoff.rd;
    [raw.cov,raw.center]=trafo(raw.cov,raw.center,med,mad,p);
    [rew.cov,rew.center]=trafo(rew.cov,rew.center,med,mad,p);
    rew.mahalanobis=abs(data'-clmean)/sqrt(clcov);
    spec.ask=1;
    %classical analysis?
    if options.classic==1
        classic.cov=clascov;
        classic.center=clasmean;
        classic.md=rew.mahalanobis;
        classic.class='COV';
    else
        classic=0;
    end
    %assigning the output
    rewo=rew;rawo=raw;
    rew=struct('center',{rewo.center},'cov',{rewo.cov},'cor',{rewo.cor},'h',{h},'Hsubsets',{rewo.Hsubsets},...
        'alpha',{alfa},'rd',{rewo.rd},'cutoff',{cutoff},'flag',{rewo.flag}, 'plane',{[]},'method',{rewo.method},...
        'class',{'MCDCOV'},'md',{rewo.mahalanobis},'classic',{classic},'X',{xx});
    raw=struct('center',{rawo.center},'cov',{rawo.cov},'cor',{rawo.cor},'objective',{rawo.objective},...
        'rd',{rawo.rd},'cutoff',{cutoff},'wt',{rawo.wt}, 'class',{'MCDCOV'},'classic',{classic},'X',{x});
    if size(data,2)~=2
        rew=rmfield(rew,'X');
        raw=rmfield(raw,'X');
    end
    try
        if plots & options.classic
            makeplot(rew,'classic',1)
        elseif plots
            makeplot(rew)
        end
    catch %output must be given even if plots are interrupted
        %> delete(gcf) to get rid of the menu
        return
    end
    return
end
%exact fit situation
if rcl < p
    % all observations lie on a hyperplane.
    z = Pcl(:,1);
    rew.plane=z;
    weights(1:n)=1;
    if cor
        correl=clcov./(sqrt(diag(clcov))*sqrt(diag(clcov))');
    end
    [clcov,clmean]=trafo(clcov,clmean,med,mad,p);
    [raw.center,raw.cov,rew.center,rew.cov,raw.objective,raw.wt,rew.flag,...
        rew.method]=displ(1,n,weights,n,p,clmean,clcov,rew.method,z./mad',ok,...
        raw.wt,rew.flag,cor,correl);
    if cor
        [rew.cor,raw.cor]=deal(correl);
    end
    %classical analysis?
    if options.classic==1
        classic.cov=clascov;
        classic.center=clasmean;
        classic.class='COV';
    else
        classic=0;
    end
    rew.Hsubsets.Hopt=1:n;
    rew.Hsubsets.Hfreq=1:n;
    %assigning the output
    rewo=rew;rawo=raw;
    rew=struct('center',{rewo.center},'cov',{rewo.cov},'cor',{rewo.cor},'h',{h},'Hsubsets',{rewo.Hsubsets},'alpha',{alfa},...
        'rd',{rewo.rd},'cutoff',{cutoff},'flag',{rewo.flag},'plane',{rewo.plane},'method',{rewo.method},...
        'class',{'MCDCOV'},'classic',{classic},'X',{xx});
    raw=struct('center',{rawo.center},'cov',{rawo.cov},'cor',{rawo.cor},'objective',{rawo.objective},...
        'cutoff',{cutoff},'wt',{rawo.wt}, 'class',{'MCDCOV'},'classic',{classic},'X',{x});
    if size(data,2)~=2
        rew=rmfield(rew,'X');
        raw=rmfield(raw,'X');
    end
    return
end

% The classical case is now handled.
if h==n
    if plots
        msg=sprintf('The MCD estimates based on %g observations are equal to the classical estimates.\n',h);
        rew.method=strvcat(rew.method,msg);
    end
    raw.center=clmean;
    raw.cov=clcov;
    raw.objective=det(clcov);
    mah=libra_mahalanobis(data,clmean,'cov',clcov);
    rew.mahalanobis=sqrt(mah);
    raw.rd=rew.mahalanobis;
    weights= mah <= cutoff.rd^2;
    raw.wt=weights;
    [rew.center,rew.cov]=weightmecov(data,weights);
    if cor
        raw.cor=raw.cov./(sqrt(diag(raw.cov))*sqrt(diag(raw.cov))');
        rew.cor=rew.cov./(sqrt(diag(rew.cov))*sqrt(diag(rew.cov))');
    else
        raw.cor=0;
        rew.cor=0;
    end
    if det(rew.cov) < exp(-50*p)
        [center,covar,z,correl,plane,count]=fit(data,NaN,med,mad,p,z,cor,rew.center,rew.cov,n);
        rew.plane=z;
        if cor
            correl=covar./(sqrt(diag(cov))*sqrt(diag(covar))');
        end
        rew.method=displrw(count,n,p,center,covar,rew.method,z,cor,correl);
        [raw.cov,raw.center]=trafo(raw.cov,raw.center,med,mad,p);
        [rew.cov,rew.center]=trafo(rew.cov,rew.center,med,mad,p);
        rew.rd=raw.rd;
    else
        mah=libra_mahalanobis(data,rew.center,'cov',rew.cov);
        weights = mah <= cutoff.md^2;
        [raw.cov,raw.center]=trafo(raw.cov,raw.center,med,mad,p);
        [rew.cov,rew.center]=trafo(rew.cov,rew.center,med,mad,p);
        rew.rd=sqrt(mah);
    end
    raw.objective=raw.objective*prod(mad)^2;
    rew.flag=weights;
    %classical analysis?
    if options.classic==1
        classic.cov=clascov;
        classic.center=clasmean;
        classic.md=rew.mahalanobis;
        classic.class='COV';
    else
        classic=0;
    end
    %assigning Hsubsets:
    rew.Hsubsets.Hopt = 1:n;
    rew.Hsubsets.Hfreq = 1:n;
    %assigning the output
    rewo=rew;rawo=raw;
    rew=struct('center',{rewo.center},'cov',{rewo.cov},'cor',{rewo.cor},'h',{h},'Hsubsets',{rewo.Hsubsets},'alpha',{alfa},...
        'rd',{rewo.rd},'cutoff',{cutoff},'flag',{rewo.flag},'plane',{rewo.plane},...
        'method',{rewo.method},'class',{'MCDCOV'},'md',{rewo.mahalanobis},'classic',{classic},'X',{xx});
    raw=struct('center',{rawo.center},'cov',{rawo.cov},'cor',{rawo.cor},'objective',{rawo.objective},...
        'rd',{rawo.rd},'cutoff',{cutoff},'wt',{rawo.wt}, 'class',{'MCDCOV'},'classic',{classic},'X',{x});
    if size(data,2)~=2
        rew=rmfield(rew,'X');
        raw=rmfield(raw,'X');
    end
    try
        if plots & options.classic
            makeplot(rew,'classic',1)
        elseif plots
            makeplot(rew)
        end
    catch %output must be given even if plots are interrupted
        %> delete(gcf) to get rid of the menu
        return
    end
    return
end
percent=h/n;
teller = zeros(1,n+1);

if Hsets_ind
    csteps = csteps1;
    inplane = NaN;
    fine = 0;
    part = 0;
    final = 1;
    tottimes = 0;
    nsamp = size(Hsets,1);
    obsingroup = n;
else

    %  If n >= 2*nmini the dataset will be divided into subdatasets.  For n < 2*nmini the set
    %  will be treated as a whole.

    if n >= 2*nmini
        maxobs=maxgroup*nmini;
        if n >= maxobs
            ngroup=maxgroup;
            group(1:maxgroup)=nmini;
        else
            ngroup=floor(n/nmini);
            minquan=floor(n/ngroup);
            group(1)=minquan;
            for s=2:ngroup
                group(s)=minquan+double(rem(n,ngroup)>=s-1);
            end
        end
        part=1;
        adjh=floor(group(1)*percent);
        nsamp=floor(ntrial/ngroup);
        minigr=sum(group);
        obsingroup=fillgroup(n,group,ngroup,seed);
        % obsingroup : i-th row contains the observations of the i-th group.
        % The last row (ngroup+1-th) contains the observations for the 2nd stage
        % of the algorithm.
    else
        [part,group,ngroup,adjh,minigr,obsingroup]=deal(0,n,1,h,n,n);
        replow=[50,22,17,15,14,zeros(1,45)];
        if n < replow(p)
            % All (p+1)-subsets will be considered.
            al=1;
            perm=[1:p,p];
            nsamp=nchoosek(n,p+1);
        else
            al=0;
            nsamp=ntrial;
        end
    end
    % some further initialisations.

    csteps=csteps1;
    inplane=NaN;
    % tottimes : the total number of iteration steps.
    % fine     : becomes 1 when the subdatasets are merged.
    % final    : becomes 1 for the final stage of the algorithm.
    [tottimes,fine,final,prevdet]=deal(0);

    if part
        % bmean1 : contains, for the first stage of the algorithm, the means of the ngroup*10
        %          best estimates.
        % bcov1  : analogous to bmean1, but now for the covariance matrices.
        % bobj1  : analogous to bmean1, but now for the objective values.
        % coeff1 : if in the k-th subdataset there are at least adjh observations that lie on
        %          a hyperplane then the coefficients of this plane will be stored in the
        %          k-th column of coeff1.
        coeff1=repmat(NaN,p,ngroup);
        bobj1=repmat(inf,ngroup,10);
        bmean1=cell(ngroup,10);
        bP1 = cell(ngroup,10);
        bL1 = cell(ngroup,10);
        [bmean1{:}]=deal(NaN);
        [bL1{:}]=deal(NaN);
        [bP1{:}]=deal(NaN);
    end

    % bmean : contains the means of the ten best estimates obtained in the second stage of the
    %         algorithm.
    % bcov  : analogous to bmean, but now for the covariance matrices.
    % bobj  : analogous to bmean, but now for the objective values.
    % coeff : analogous to coeff1, but now for the merged subdataset.
    % If the data is not split up, the 10 best estimates obtained after csteps1 iterations
    % will be stored in bmean, bcov and bobj.
    coeff=repmat(NaN,p,1);
    bobj=repmat(inf,1,10);
    bmean=cell(1,10);
    bP = cell(1,10);
    bL = cell(1,10);
    [bmean{:}]=deal(NaN);
    [bP{:}]=deal(NaN);
    [bL{:}] = deal(NaN);
end

seed=0;
tellerwhilelus = 0;

while final~=2
    if fine | (~part & final)
        if ~Hsets_ind
            nsamp=10;
        end
        if final
            adjh=h;
            ngroup=1;
            if n*p <= 1e+5
                csteps=csteps3;
            elseif n*p <=1e+6
                csteps=10-(ceil(n*p/1e+5)-2);
            else
                csteps=1;
            end
            if n > 5000

                nsamp=1;
            end
        else
            adjh=floor(minigr*percent);
            csteps=csteps2;
        end
    end

    % found : becomes 1 if we have a singular intermediate MCD estimate.
    found=0;

    for k=1:ngroup
        if ~fine
            found=0;
        end
        for i=1:nsamp
            tottimes=tottimes+1;
            % ns becomes 1 if we have a singular trial subsample and if there are at
            % least adjh observations in the subdataset that lie on the concerning hyperplane.
            % In that case we don't have to take C-steps. The determinant is zero which is
            % already the lowest possible value. If ns=1, no C-steps will be taken and we
            % start with the next sample. If we, for the considered subdataset, haven't
            % already found a singular MCD estimate, then the results must be first stored in
            % bmean, bcov, bobj or in bmean1, bcov1 and bobj1.  If we, however, already found
            % a singular result for that subdataset, then the results won't be stored
            % (the hyperplane we just found is probably the same as the one we found earlier.
            % We then let adj be zero. This will guarantee us that the results won't be
            % stored) and we start immediately with the next sample.
            adj=1;
            ns=0;

            % For the second and final stage of the algorithm the array sortdist(1:adjh)
            % contains the indices of the observations corresponding to the adjh observations
            % with minimal relative distances with respect to the best estimates of the
            % previous stage. An exception to this, is when the estimate of the previous
            % stage is singular.  For the second stage we then distinguish two cases :
            %
            % 1. There aren't adjh observations in the merged set that lie on the hyperplane.
            %    The observations on the hyperplane are then extended to adjh observations by
            %    adding the observations of the merged set with smallest orthogonal distances
            %    to that hyperplane.
            % 2. There are adjh or more observations in the merged set that lie on the
            %    hyperplane. We distinguish two cases. We haven't or have already found such
            %    a hyperplane. In the first case we start with a new sample.  But first, we
            %    store the results in bmean1, bcov1 and bobj1. In the second case we
            %    immediately start with a new sample.
            %
            % For the final stage we do the same as 1. above (if we had h or more observations
            % on the hyperplane we would already have found it).

            if ~Hsets_ind

                if final
                    if ~isinf(bobj(i))
                        meanvct=bmean{i};
                        P=bP{i};
                        L = bL{i};
                        if bobj(i)==0
                            [dis,sortdist]=sort(abs(sum((data-repmat(meanvct,n,1))'.*repmat(coeff,1,n))));
                        else
                            [dis,sortdist]=mahal2((data-repmat(meanvct,n,1))*P,sqrt(L),part,fine,final,k,obsingroup);
                        end
                    else
                        break
                    end
                elseif fine
                    if ~isinf(bobj1(k,i))
                        meanvct=bmean1{k,i};
                        P=bP1{k,i};
                        L=bL1{k,i};
                        if bobj1(k,i)==0
                            [dis,ind]=sort(abs(sum((data(obsingroup{end},:)-repmat(meanvct,minigr,1))'.*repmat(coeff1(:,k),1,minigr))));
                            sortdist=obsingroup{end}(ind);
                            if dis(adjh) < 1e-8
                                if found==0
                                    obj=0;
                                    coeff=coeff1(:,k);
                                    found=1;
                                else
                                    adj=0;
                                end
                                ns=1;
                            end
                        else
                            [dis,sortdist]=mahal2((data(obsingroup{end},:)-repmat(meanvct,minigr,1))*P,sqrt(L),part,fine,final,k,obsingroup);
                        end
                    else
                        break;
                    end
                else
                    if ~part
                        if al
                            k=p+1;
                            perm(k)=perm(k)+1;
                            while ~(k==1 |perm(k) <=(n-(p+1-k)))
                                k=k-1;
                                perm(k)=perm(k)+1;
                                for j=(k+1):p+1
                                    perm(j)=perm(j-1)+1;
                                end
                            end
                            index=perm;
                        else
                            [index,seed]=randomset(n,p+1,seed);
                        end
                    else
                        [index,seed]=randomset(group(k),p+1,seed);
                        index=obsingroup{k}(index);
                    end

                    [P,T,L,r,centerX,meanvct] = classSVD(data(index,:));
                    if r==0
                        ns = 1;
                        rew.center=meanvct;
                        rew.cov=zeros(p,p);
                        rew.flag = zeros(1,n);
                        rew.flag(index) = 1;
                        rew.Hsubsets.Hopt=1:n;
                        rew.Hsubsets.Hfreq=1:n;
                    elseif (r < p) & (r~=0)
                        % The trial subsample is singular.
                        % We distinguish two cases :
                        %
                        % 1. There are adjh or more observations in the subdataset that lie
                        %    on the hyperplane. If the data is not split up, we have adjh=h and thus
                        %    an exact fit. If the data is split up we distinguish two cases.
                        %    We haven't or have already found such a hyperplane.  In the first case
                        %    we check if there are more than h observations in the entire set
                        %    that lie on the hyperplane. If so, we have an exact fit situation.
                        %    If not, we start with a new trial subsample.  But first, the
                        %    results must be stored bmean1, bcov1 and bobj1.  In the second case
                        %    we immediately start with a new trial subsample.
                        %
                        % 2. There aren't adjh observations in the subdataset that lie on the
                        %    hyperplane. We then extend the trial subsample until it isn't singular
                        %    anymore.


                        % eigvct : contains the coefficients of the hyperplane.

                        eigvct = P(:,1);

                        if ~part
                            dist=abs(sum((data-repmat(meanvct,n,1))'.*repmat(eigvct,1,n)));
                        else
                            dist=abs(sum((data(obsingroup{k},:)-repmat(meanvct,group(k),1))'.*repmat(eigvct,1,group(k))));
                        end

                        obsinplane=find(dist < 1e-8);
                        % count : number of observations that lie on the hyperplane.
                        count=length(obsinplane);

                        if count >= adjh
                            if ~part
                                [center,covar,eigvct,correl]=fit(data,obsinplane,med,mad,p,eigvct,cor);
                                rew.plane=eigvct;
                                weights(obsinplane)=1;
                                [raw.center,raw.cov,rew.center,rew.cov,raw.objective,...
                                    raw.wt,rew.flag,rew.method]=displ(2,count,weights,n,p,center,covar,...
                                    rew.method,eigvct,ok,raw.wt,rew.flag,cor,correl);
                                if cor
                                    [rew.cor,raw.cor]=deal(correl);
                                end
                                rew.Hsubsets.Hopt=obsinplane;
                                rew.Hsubsets.Hfreq=obsinplane;
                                return
                            elseif found==0
                                if count>=h
                                    [center,covar,eigvct,correl]=fit(data,obsinplane,med,mad,p,eigvct,cor);
                                    rew.plane=eigvct;
                                    weights(obsinplane)=1;
                                    [raw.center,raw.cov,rew.center,rew.cov,raw.objective,...
                                        raw.wt,rew.flag,rew.method,varargout]=displ(2,count2,weights,n,p,center,covar,...
                                        rew.method,eigvct,ok,raw.wt,rew.flag,cor,correl);
                                    if cor
                                        [rew.cor,raw.cor]=deal(correl);
                                    end
                                    rew.Hsubsets.Hopt=obsinplane;
                                    rew.Hsubsets.Hfreq=obsinplane;
                                    return
                                end
                                obj=0;
                                inplane(k)=count;
                                coeff1(:,k)=eigvct;
                                found=1;
                                ns=1;
                            else
                                ns=1;
                                adj=0;
                            end
                        else
                            covmat = cov(data(index,:));
                            meanvct = mean(data(index,:));
                            while det(covmat) < exp(-50*p)
                                [index1,seed]=addobs(index,n,seed);
                                [covmat,meanvct] = updatecov(data(index,:),covmat,meanvct,data(setdiff(index1,index),:),[],1);
                                index = index1;
                            end
                        end
                    end

                    if ~ns
                        if ~part
                            [dis,sortdist] = mahal2((data - repmat(meanvct,n,1))*P,sqrt(L),part,fine,final,k,obsingroup);
                        else
                            [dis,sortdist] = mahal2((data(obsingroup{k},:) - repmat(meanvct,group(k),1))*P,sqrt(L),part,fine,final,k,obsingroup);
                        end
                    end
                end
            end

            if ~ns
                for j=1:csteps
                    tottimes=tottimes+1;
                    if j == 1
                        if Hsets_ind
                            obs_in_set = Hsets(i,:);
                        else
                            obs_in_set = sort(sortdist(1:adjh));
                            teller(obs_in_set) = teller(obs_in_set) + 1;
                            teller(end) = teller(end) + 1;
                        end
                    else
                        % The observations correponding to the adjh smallest mahalanobis
                        % distances determine the subset for the next iteration.
                        if ~part
                            [dis2,sortdist] = mahal2((data - repmat(meanvct,n,1))*P,sqrt(L),part,fine,final,k,obsingroup);
                        else
                            if final
                                [dis2,sortdist] = mahal2((data - repmat(meanvct,n,1))*P,sqrt(L),part,fine,final,k,obsingroup);
                            elseif fine
                                [dis2,sortdist] = mahal2((data(obsingroup{end},:) - repmat(meanvct,minigr,1))*P,sqrt(L),part,fine,final,k,obsingroup);
                            else
                                [dis2,sortdist] = mahal2((data(obsingroup{k},:) - repmat(meanvct,group(k),1))*P,sqrt(L),part,fine,final,k,obsingroup);
                            end
                        end
                        % Creation of a H-subset.
                        obs_in_set=sort(sortdist(1:adjh));
                        teller(obs_in_set) = teller(obs_in_set) + 1;
                        teller(end) = teller(end) + 1;
                    end
                    [P,T,L,r,centerX,meanvct] = classSVD(data(obs_in_set,:));
                    if r==0
                        rew.center=meanvct;
                        rew.cov=zeros(p,p);
                        rew.flag=(data==data(obs_in_set(1),:));
                        rew.Hsubset.Hopt=obs_in_set(1);
                        rew.Hsubset.Hfreq=obs_in_set(1);
                    else

                        obj=prod(L);

                        if obj < exp(-50*p)
                            % The adjh-subset is singular. If adjh=h we have an exact fit situation.
                            % If adjh < h we distinguish two cases :
                            %
                            % 1. We haven't found earlier a singular adjh-subset. We first check if
                            %    in the entire set there are h observations that lie on the hyperplane.
                            %    If so, we have an exact fit situation. If not, we stop taking C-steps
                            %    (the determinant is zero which is the lowest possible value) and
                            %    store the results in the appropriate arrays.  We then begin with
                            %    the next trial subsample.
                            %
                            % 2. We have, for the concerning subdataset, already found a singular
                            %    adjh-subset. We then immediately begin with the next trial subsample.

                            if ~part | final | (fine & n==minigr)
                                covmat=cov(data(obs_in_set,:));
                                [center,covar,z,correl,obsinplane,count]=fit(data,NaN,med,mad,p,NaN,...
                                    cor,meanvct,covmat,n);
                                rew.plane=z;
                                weights(obsinplane)=1;
                                [raw.center,raw.cov,rew.center,rew.cov,raw.objective,...
                                    raw.wt,rew.flag,rew.method]=displ(2,count,weights,n,p,center,covar,...
                                    rew.method,z,ok,raw.wt,rew.flag,cor,correl);
                                if cor
                                    [rew.cor,raw.cor]=deal(correl);
                                end
                                rew.Hsubset.Hopt=obsinplane;
                                rew.Hsubset.Hfreq=obsinplane;
                                return
                            elseif found==0
                                eigvct = V(:,1);
                                dist=abs(sum((data-repmat(meanvct,n,1))'.*repmat(eigvct,1,n)));
                                obsinplane=find(dist<1e-8);
                                count=length(obsinplane);
                                if count >= h
                                    [center,covar,eigvct,correl]=fit(data,obsinplane,med,mad,p,eigvct,cor);
                                    rew.plane=eigvct;
                                    weights(obsinplane)=1;
                                    [raw.center,raw.cov,rew.center,rew.cov,raw.objective,...
                                        raw.wt,rew.flag,rew.method]=displ(2,count,weights,n,p,center,covar,...
                                        rew.method,eigvct,ok,raw.wt,rew.flag,cor,correl);
                                    if cor
                                        [rew.cor,raw.cor]=deal(correl);
                                    end
                                    rew.Hsubset.Hopt=obsinplane;
                                    rew.Hsubset.Hfreq=obsinplane;
                                    return
                                end
                                obj=0;
                                found=1;
                                if ~fine
                                    coeff1(:,k)=eigvct;
                                    dist=abs(sum((data(obsingroup{k},:)-repmat(meanvct,group(k),1))'.*repmat(eigvct,1,group(k))));
                                    inplane(k)=length(dist(dist<1e-8));
                                else
                                    coeff=eigvct;
                                    dist=abs(sum((data(obsingroup{end},:)-repmat(meanvct,minigr,1))'.*repmat(eigvct,1,minigr)));
                                    inplane=length(dist(dist<1e-8));
                                end
                                break;
                            else
                                adj=0;
                                break;
                            end
                        end
                    end
                    % We stop taking C-steps when two subsequent determinants become equal.
                    % We have then reached convergence.
                    if j >= 2 & obj == prevdet
                        break;
                    end
                    prevdet=obj;

                end % C-steps

            end

            % After each iteration, it has to be checked whether the new solution
            % is better than some previous one.  A distinction is made between the
            % different stages of the algorithm:
            %
            %  - Let us first consider the first (second) stage of the algorithm.
            %    We distinguish two cases if the objective value is lower than the largest
            %    value in bobj1 (bobj) :
            %
            %      1. The new objective value did not yet occur in bobj1 (bobj).  We then store
            %         this value, the corresponding mean and covariance matrix at the right
            %         place in resp. bobj1 (bobj), bmean1 (bmean) and bcov1 (bcov).
            %         The objective value is inserted by shifting the greater determinants
            %         upwards. We perform the same shifting in bmean1 (bmean) and bcov1 (bcov).
            %
            %      2. The new objective value already occurs in bobj1 (bobj). A comparison is
            %         made between the new mean vector and covariance matrix and those
            %         estimates with the same determinant. When for an equal determinant,
            %         the mean vector or covariance matrix do not correspond, the new results
            %         will be stored in bobj1 (bobj), bmean1 (bmean) and bcov1 (bcov).
            %
            %    If the objective value is not lower than the largest value in bobj1 (bobj),
            %    nothing happens.
            %
            %  - For the final stage of the algorithm, only the best solution has to be kept.
            %    We then check if the objective value is lower than the till then lowest value.
            %    If so, we have a new best solution. If not, nothing happens.


            if ~final & adj
                if fine | ~part
                    if obj < max(bobj) & ~ns
                        [bmean,bP,bL,bobj]=insertion(bmean,bP,bL,bobj,meanvct,P,L,obj,1,eps);
                    end
                else
                    if obj < max(bobj1(k,:)) & ~ns
                        [bmean1,bP1,bL1,bobj1]=insertion(bmean1,bP1,bL1,bobj1,meanvct,P,L,obj,k,eps);
                    end
                end
            end

            if final & obj< bestobj
                % bestset           : the best subset for the whole data.
                % bestobj           : objective value for this set.
                % initmean, initcov : resp. the mean and covariance matrix of this set.
                bestset=obs_in_set;
                bestobj=obj;
                initmean=meanvct;
                initcov=cov(data(bestset,:));
                raw.initcov=cov(data(bestset,:));
            end

        end % nsamp
    end % ngroup


    if part & ~fine
        fine=1;
    elseif (part & fine & ~final) | (~part & ~final)
        final=1;
    else
        final=2;
    end

end % while loop

[P,T,L,r,centerX,cX] = classSVD(data(bestset,:));
mah=libra_mahalanobis((data - repmat(cX,n,1))*P,zeros(size(P,2),1),'cov',L);
sortmah=sort(mah);

[sortset,indbestset] = sort(mah(bestset));
sortbestset = bestset(indbestset);
rew.Hsubsets.Hopt = sortbestset;

if ~factor_ind
    factor = sortmah(h)/chi2inv(h/n,p);
else
    factor = sortmah(h)/chi2inv(h/n,p/2);
end

raw.cov=factor*initcov;
% We express the results in the original units.
[raw.cov,raw.center]=trafo(raw.cov,initmean,med,mad,p);
raw.objective=bestobj*prod(mad)^2;

if cor
    raw.cor=raw.cov./(sqrt(diag(raw.cov))*sqrt(diag(raw.cov))');
end

% the mahalanobis distances are computed without the factor, therefore we
% have to correct for it now.
mah=mah/factor;
raw.rd=sqrt(mah);
weights=mah<=cutoff.md^2;
raw.wt=weights;
[rew.center,rew.cov]=weightmecov(data,weights);
[trcov,trcenter]=trafo(rew.cov,rew.center,med,mad,p);

% determination of Hfreq:
[telobs,indobs] = greatsort(teller(1:(end - 1)));
rew.Hsubsets.Hfreq = indobs(1:(h));
if size(rew.Hsubsets.Hfreq,2) == 1
    rew.Hsubsets.Hfreq = rew.Hsubsets.Hfreq';
end

if cor
    rew.cor=rew.cov./(sqrt(diag(rew.cov))*sqrt(diag(rew.cov))');
end

if prod(sqrt(L)) < exp(-50*p)
    [center,covar,z,correl,plane,count]=fit(data,NaN,med,mad,p,z,cor,rew.center,rew.cov,n);
    rew.plane=z;
    if cor
        correl=covar./(sqrt(diag(covar))*sqrt(diag(covar))');
        [rew.cor,raw.cor] = deal(correl);
    end
    rew.method=displrw(count,n,p,center,covar,rew.method,z,cor,correl);
    rew.flag=weights;
    rew.rd=raw.rd;
else
    mah=libra_mahalanobis(data,rew.center,'cov',rew.cov);
    rew.flag=(mah <= cutoff.md^2);
    rew.rd=sqrt(mah);
end

rew.mahalanobis=sqrt(libra_mahalanobis(data,clmean,'cov',clcov));
rawo=raw;
reso=rew;
if options.classic==1
    classic.cov=clascov;
    classic.center=clasmean;
    classic.md=rew.mahalanobis;
    classic.flag = (classic.md <= cutoff.md);
    if options.cor==1
        classic.cor=clascov./(sqrt(diag(clascov))*sqrt(diag(clascov))');
    end
    classic.class='COV';
else
    classic=0;
end
%assigning the output
rew=struct('center',{trcenter},'cov',{trcov},'cor',{reso.cor},'h',{h},'Hsubsets',{reso.Hsubsets},'alpha',{alfa},...
    'rd',{reso.rd},'flag',{reso.flag},'md',{reso.mahalanobis},'cutoff',{cutoff},...
    'plane',{reso.plane},'method',{reso.method},'class',{'MCDCOV'},'classic',{classic},'X',{xx});
raw=struct('center',{rawo.center},'cov',{rawo.cov},'cor',{rawo.cor},'objective',{rawo.objective},...
    'rd',{rawo.rd},'cutoff',{cutoff},...
    'wt',{rawo.wt},'class',{'MCDCOV'},'classic',{classic},'X',{x});

if size(data,2)~=2
    rew=rmfield(rew,'X');
    raw=rmfield(raw,'X');
end

try
    if plots & options.classic
        makeplot(rew,'classic',1)
    elseif plots
        makeplot(rew)
    end
catch %output must be given even if plots are interrupted
    %> delete(gcf) to get rid of the menu
end

%-----------------------------------------------------------------------------------------
function [raw_center,raw_cov,center,covar,raw_objective,raw_wt,mcd_wt,method]=displ(exactfit,...
    count,weights,n,p,center,covar,method,z,ok,raw_wt,mcd_wt,cor,correl,varargin)

% Determines some fields of the output argument REW for the exact fit situation.  It also
% displays and writes the messages concerning the exact fit situation.  If the raw MCD
% covariance matrix is not singular but the reweighted is, then the function displrw is
% called instead of this function.

[raw_center,center]=deal(center);
[raw_cov,cover]=deal(covar);
raw_objective=0;
mcd_wt=weights;
raw_wt=weights;

switch exactfit
    case 1
        msg='The covariance matrix of the data is singular.';
    case 2
        msg='The covariance matrix has become singular during the iterations of the MCD algorithm.';
    case 3
        msg=sprintf('The %g-th order statistic of the absolute deviation of variable %g is zero. ',varargin{1},varargin{2});
end

msg=sprintf([msg '\nThere are %g observations in the entire dataset of %g observations that lie on the \n'],count,n);
switch p
    case 2
        msg=sprintf([msg 'line with equation %g(x_i1-m_1)%+g(x_i2-m_2)=0 \n'],z);
        msg=sprintf([msg 'where the mean (m_1,m_2) of these observations is the MCD location']);
    case 3
        msg=sprintf([msg 'plane with equation %g(x_i1-m_1)%+g(x_i2-m_2)%+g(x_i3-m_3)=0 \n'],z);
        msg=sprintf([msg 'where the mean (m_1,m_2,m_3) of these observations is the MCD location']);
    otherwise
        msg=sprintf([msg 'hyperplane with equation a_1 (x_i1-m_1) + ... + a_p (x_ip-m_p) = 0 \n']);
        msg=sprintf([msg 'with coefficients a_i equal to : \n\n']);
        msg=sprintf([msg sprintf('%g  ',z)]);
        msg=sprintf([msg '\n\nand where the mean (m_1,...,m_p) of these observations is the MCD location']);
end

method=strvcat(method,[msg '.']);
disp(method)

%-----------------------------------------------------------------------------------------
function method=displrw(count,n,p,center,covar,method,z,cor,correl)

% Displays and writes messages in the case the reweighted robust covariance matrix
% is singular.

msg=sprintf('The reweighted MCD scatter matrix is singular. \n');
msg=sprintf([msg 'There are %g observations in the entire dataset of %g observations that lie on the\n'],count,n);

switch p
    case 2
        msg=sprintf([msg 'line with equation %g(x_i1-m_1)%+g(x_i2-m_2)=0 \n\n'],z);
        msg=sprintf([msg 'where the mean (m_1,m_2) of these observations is : \n\n']);
    case 3
        msg=sprintf([msg 'plane with equation %g(x_i1-m_1)%+g(x_i2-m_2)%+g(x_i3-m_3)=0 \n\n'],z);
        msg=sprintf([msg 'where the mean (m_1,m_2,m_3) of these observations is : \n\n']);
    otherwise
        msg=sprintf([msg 'hyperplane with equation a_1 (x_i1-m_1) + ... + a_p (x_ip-m_p) = 0 \n']);
        msg=sprintf([msg 'with coefficients a_i equal to : \n\n']);
        msg=sprintf([msg sprintf('%g  ',z)]);
        msg=sprintf([msg '\n\nand where the mean (m_1,...,m_p) of these observations is : \n\n']);
end

msg=sprintf([msg sprintf('%g  ',center)]);
msg=sprintf([msg '\n\nTheir covariance matrix equals : \n\n']);
msg=sprintf([msg sprintf([repmat('% 13.5g ',1,p) '\n'],covar)]);
if cor
    msg=sprintf([msg '\n\nand their correlation matrix equals : \n\n']);
    msg=sprintf([msg sprintf([repmat('% 13.5g ',1,p) '\n'],correl)]);
end

method=strvcat(method,msg);

%-----------------------------------------------------------------------------------------
function [initmean,initcov,z,correl,varargout]=fit(dat,plane,med,mad,p,z,cor,varargin)

% This function is called in the case of an exact fit. It computes the correlation
% matrix and transforms the coefficients of the hyperplane, the mean, the covariance
% and the correlation matrix to the original units.

if isnan(plane)
    [meanvct,covmat,n]=deal(varargin{:});
    [z, eigvl]=eigs(covmat,1,0,struct('disp',0));
    dist=abs(sum((dat-repmat(meanvct,n,1))'.*repmat(z,1,n)));
    plane=find(dist < 1e-8);
    varargout{1}=plane;
    varargout{2}=length(plane);
end

z=z./mad';
[initcov,initmean]=trafo(cov(dat(plane,:)),mean(dat(plane,:)),med,mad,p);
if cor
    correl=initcov./(sqrt(diag(initcov))*sqrt(diag(initcov))');
else
    correl=NaN;
end

%------------------------------------------------------------------------------------------
function obsingroup=fillgroup(n,group,ngroup,seed)

% Creates the subdatasets.

obsingroup=cell(1,ngroup+1);

jndex=0;
for k=1:ngroup
    for m=1:group(k)
        [random,seed]=uniran(seed);
        ran=floor(random*(n-jndex)+1);
        jndex=jndex+1;
        if jndex==1
            index(1,jndex)=ran;
            index(2,jndex)=k;
        else
            index(1,jndex)=ran+jndex-1;
            index(2,jndex)=k;
            ii=min(find(index(1,1:jndex-1) > ran-1+[1:jndex-1]));
            if length(ii)
                index(1,jndex:-1:ii+1)=index(1,jndex-1:-1:ii);
                index(2,jndex:-1:ii+1)=index(2,jndex-1:-1:ii);
                index(1,ii)=ran+ii-1;
                index(2,ii)=k;
            end
        end
    end
    obsingroup{k}=index(1,index(2,:)==k);
    obsingroup{ngroup+1}=[obsingroup{ngroup+1},obsingroup{k}];
end

%-----------------------------------------------------------------------------------------
function [ranset,seed]=randomset(tot,nel,seed)

% This function is called if not all (p+1)-subsets out of n will be considered.
% It randomly draws a subsample of nel cases out of tot.

for j=1:nel
    [random,seed]=uniran(seed);
    num=floor(random*tot)+1;
    if j > 1
        while any(ranset==num)
            [random,seed]=uniran(seed);
            num=floor(random*tot)+1;
        end
    end
    ranset(j)=num;
end

%-----------------------------------------------------------------------------------------
function [index,seed]=addobs(index,n,seed)

% Extends a trial subsample with one observation.

jndex=length(index);
[random,seed]=uniran(seed);
ran=floor(random*(n-jndex)+1);
jndex=jndex+1;
index(jndex)=ran+jndex-1;
ii=min(find(index(1:jndex-1) > ran-1+[1:jndex-1]));
if length(ii)~=0
    index(jndex:-1:ii+1)=index(jndex-1:-1:ii);
    index(ii)=ran+ii-1;
end

%-----------------------------------------------------------------------------------------

function mahsort=mahal(dat,meanvct,covmat,part,fine,final,k,obsingroup,group,minigr,n,nvar)

% Orders the observations according to the mahalanobis distances.

if ~part | final
    [dis,ind]=sort(libra_mahalanobis(dat,meanvct,'cov',covmat));
    mahsort=ind;
elseif fine
    [dis,ind]=sort(libra_mahalanobis(dat(obsingroup{end},:),meanvct,'cov',covmat));
    mahsort=obsingroup{end}(ind);
else
    [dis,ind]=sort(libra_mahalanobis(dat(obsingroup{k},:),meanvct,'cov',covmat));
    mahsort=obsingroup{k}(ind);
end

%-----------------------------------------------------------------------------------------

function [dis,mahsort]=mahal2(score,sca,part,fine,final,k,obsingroup)

% Orders the observations according to the mahalanobis distances for a diagonal
% covariance matrix and zero mean. sca contains the squareroot of the diagonal elements.

if ~part | final
    [dis,ind]=sort(libra_mahalanobis(score,zeros(size(score,2),1),'cov',sca.^2));
    mahsort=ind;
elseif fine
    [dis,ind]=sort(libra_mahalanobis(score,zeros(size(score,2),1),'cov',sca.^2));
    mahsort=obsingroup{end}(ind);
else
    [dis,ind]=sort(libra_mahalanobis(score,zeros(size(score,2),1),'cov',sca.^2));
    mahsort=obsingroup{k}(ind);
end

%------------------------------------------------------------------------------------------

function [covmat,meanvct]=trafo(covmat,meanvct,med,mad,nvar)

% Transforms a mean vector and a covariance matrix to the original units.

covmat=covmat.*repmat(mad,nvar,1).*repmat(mad',1,nvar);
meanvct=meanvct.*mad+med;

%-----------------------------------------------------------------------------------------
function [bestmean,bestP,bestL,bobj]=insertion(bestmean,bestP,bestL,bobj,meanvct,P,L,obj,row,eps)

% Stores, for the first and second stage of the algorithm, the results in the appropriate
% arrays if it belongs to the 10 best results.

insert=1;

equ=find(obj==bobj(row,:));

for j=equ
    if (meanvct==bestmean{row,j})
        if all(P==bestP{row,j})
            if all(L==bestL{row,j})
                insert=0;
            end
        end
    end
end

if insert
    ins=min(find(obj < bobj(row,:)));

    if ins==10
        bestmean{row,ins}=meanvct;
        bestP{row,ins} = P;
        bestL{row,ins} = L;
        bobj(row,ins)=obj;
    else
        [bestmean{row,ins+1:10}]=deal(bestmean{row,ins:9});
        bestmean{row,ins}=meanvct;
        [bestP{row,ins+1:10}] = deal(bestP{row,ins:9});
        bestP{row,ins} = P;
        [bestL{row,ins+1:10}] = deal(bestL{row,ins:9});
        bestL{row,ins} = L;
        bobj(row,ins+1:10)=bobj(row,ins:9);
        bobj(row,ins)=obj;
    end

end

%-----------------------------------------------------------------------------------------
function quan=quanf(alfa,n,rk)
quan=floor(2*floor((n+rk+1)/2)-n+2*(n-floor((n+rk+1)/2))*alfa);
%--------------------------------------------------------------------------
