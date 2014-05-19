function [rew,raw] = ltsregres(x,y,varargin);

%LTSREGRES carries out least trimmed squares (LTS) regression, introduced in
%   
%   Rousseeuw, P.J. (1984), "Least Median of Squares Regression,"
%   Journal of the American Statistical Association, Vol. 79, pp. 871-881.
% 
% The LTS regression method minimizes the sum of the h smallest squared 
% residuals, where h must be at least half the number of observations. The 
% default value of h is roughly 0.75n (n is the total number of observations), 
% but the user may choose any value between n/2 and n.
%
% To compute the LTS estimator, the FAST-LTS algorithm is used.
% Reference: 
%  Rousseeuw, P.J. and Van Driessen, K. (2006), 
%  "Computing LTS Regression for Large Data Sets", Data Mining and
%  Knowledge Discovery, 12, 29-45.
%
%  Also available at http://www.agoras.ua.ac.be/
%
% The LTS regression method is intended for continuous variables, and assumes 
% that the number of observations n is at least 2 times the number of 
% regression coefficients p. If p is too large with respect to n, it is  
% better to first reduce p by variable selection or principal components 
% (see rpcr.m, rsimpls.m). The response variable should be univariate, otherwise 
% robust multivariate regression should be performed (see mcdregres.m).
%
% The LTS is a robust method in the sense that the estimated regression
% fit is not unduly influenced by outliers in the data, even if there are
% several outliers. Due to this robustness, we can detect outliers by their
% large LTS residuals.
%
% Required input arguments: 
%    x : Data matrix of explanatory variables (also called 'regressors'). 
%        Rows of x represent observations, and columns represent variables.  
%        Missing values (NaN's) and infinite values (Inf's) are allowed, since observations (rows) 
%        with missing or infinite values will automatically be excluded from the computations.
%    y:  A vector with n elements that contains the response variables.
%        Missing values (NaN's) and infinite values (Inf's) are allowed, since observations (rows) 
%        with missing or infinite values will automatically be excluded from the computations.
%
% Optional input arguments:           
%   intercept : If 1, a model with constant term will be fitted (default),
%               if 0, no constant term will be included. 
%   intadjust : If 1, the intercept adjustment will be applied in each step
%               of the algorithm. These calculations need substantially more 
%               computation time than intadjust=0, which is the default value.
%           h : The number of observations that have determined the least 
%               trimmed squares estimator. Any value between n/2 and n may be specified.
%       alpha : (1-alpha) measures the fraction of outliers the algorithm should 
%               resist. Any value between 0.5 and 1 may be specified. (default = 0.75)
%      ntrial : Number of initial subsets drawn. Its default value is 500.
%       plots : If equal to one, a menu is shown which allows to draw several plots,
%               such as residual plots and a regression outlier map. (default)
%               If the input argument 'classic' is equal to one, the classical
%               plots are drawn as well.
%               If 'plots' is equal to zero, all plots are suppressed.
%               See also makeplot.m
%     classic : If equal to one, classical least squares regression will be performed,
%               see ols.m (default = 0).
%
% Input arguments for advanced users:
%     Hsets : Instead of random trial h-subsets (default, Hsets = []), Hsets makes it possible to give certain
%             h-subsets as input. Hsets is a matrix that contains the indices of the observations of one
%             h-subset as a row.
%
% I/O:  result=ltsregres(x,y,'plots',0,'intercept',0)
%       [rew,raw] = ltsregres(x,y)
%  The user should only give the input arguments that have to change their default value.
%  The name of the input arguments needs to be followed by their value.
%  The order of the input arguments is of no importance.
%
% The output consists of two structures 'rew' and 'raw' containing the 
% following fields:
%
%    raw.coefficients : Vector of raw LTS coefficient estimates (including the 
%                       intercept, when options.intercept=1).
%          raw.fitted : Vector like y containing the raw fitted values of the response.
%             raw.res : Vector like y containing the raw residuals from the regression.
%           raw.scale : Scale estimate of the raw residuals.
%       raw.objective : Objective function of the LTS regression method, i.e. the sum 
%                       of the h smallest squared raw residuals.
%              raw.wt : Vector like y containing weights that can be used in a weighted
%                       least squares. These weights are 1 for points with reasonably
%                       small raw residuals, and 0 for points with large raw residuals.
%
% 	        rew.slope : Vector of the slope coefficients obtained after reweighting.
%             rew.int : The intercept.
%          rew.fitted : Vector like y containing the fitted values of the response
%                       after reweighting.
%             rew.res : Vector like y containing the residuals from the weighted
%                       least squares regression.
%           rew.scale : Scale estimate of the reweighted residuals.
%        rew.rsquared : Robust version of R squared. This is 1 minus the fraction:
%                       (sum of the quan smallest squared residuals) over (sum of 
%                       the quan smallest (y-loc)^2), where the denominator
%                       is minimized over loc. Note that loc is not subtracted from
%                       y if intercept = 0 in the call to ltsregres.	
%               rew.h : The number of observations that have determined the LTS estimator, 
%                       i.e. the value of h. 
%       rew. Hsubsets : A structure that contains Hopt and Hfreq:
%                        Hopt  : The subset of h points whose covariance matrix has minimal determinant, 
%                                ordered following increasing robust distances.
%                        Hfreq : The subset of h points which are the most frequently selected during the whole
%                                algorithm.
%           rew.alpha : (1-alpha) measures the fraction of outliers the algorithm should 
%                       resist. 
%              rew.rd : The robust distances for the observations of the design matrix, 
%                       based on the MCD estimator (mcdcov.m)  
%            rew.resd : Vector like y containing the standardized residuals
%                       from the weighted least squares regression.
%         rew.weights : Vector like y containing weights that have been used in a weighted
%                       least squares. These weights are 1 for points with reasonably
%                       small raw residuals, and 0 for points with large raw residuals.
%          rew.cutoff : Structure which contains cutoff values for the robust distances computed by mcdcov.m,
%                       and for the standardized residuals. 
%            rew.flag : Vector like y containing flags based on the reweighted regression.
%                       These flags determine which observations can be considered as 
%                       outliers.
%          rew.method : Character string naming the method (Least Trimmed Squares).	
%           rew.class : 'LTS'
%         rew.classic : If the input argument 'classic' equals 1, this structure contains the 
%                       results of a classical least squares regression
%               rew.X : If x is univariate, same as the input x in the call to ltsregres,
%                       without rows containing missing or infinite values. 
%               rew.y : If x is univariate, same as the input y in the call to ltsregres,
%                       without rows containing missing or infinite values.

%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Version 22/12/2000, 
% Written by Katrien Van Driessen and Randy Brenkers
% Revisions by Sabine Verboven, Sanne Engelen, Nele Smets
% Last update: 06/07/2004

%MATLAB: to avoid singularMatrix in line 538 
warning off 
%The maximum value for p = number of variables
pmax=20;
%The maximum value for n = number of observations
nmax=50000;
%To change the number of subdatasets and their size, the values of maxgroup
%and nmini can be changed
maxgroup=5;
nmini=300;
%The number of iteration steps in stages 1,2 and 3 can be changed
% by adapting the parameters csteps1, csteps2, and csteps3.
csteps1=2;
csteps2=2;
csteps3=100;

pmax1=pmax+1;
pmax2=pmax*pmax;
nvm11=pmax*pmax1;
nvm12=pmax1*pmax1;
km10=10*maxgroup;
nmaxi=nmini*maxgroup;
maxmini=fix(((3*nmini-1)/2)+1);
% dtrial : number of subsamples if not all (p+1)-subsets will be considered. 
dtrial=500;

[n,p]=size(x);
[m,q]=size(y);
if q~=1
    if m==1
        y=y';
    else
        error('y is not one-dimensional.');
    end
end
na.x=~isfinite(x*ones(p,1));
na.y=~isfinite(y);
if size(na.x,1)~=size(na.y,1)
    error('Number of observations in x and y not equal.');
end
% Observations with missing or infinite values are ommitted. 
ok=~(na.x|na.y);
x=x(ok,:);
y=y(ok,:);
dx=size(x);
dy=size(y,1);
n=dx(1);
% Some checks are now performed.
if n == 0
    error('All observations have missing values!');
end
if n > nmax
    error(['The program allows for at most ' int2str(nmax) ' observations.']);
end
%internal variables and default values
seed=0;
alfa=0.75;
hdef=quanf(alfa,n,p+1); %default value of h
hmin=quanf(0.5,n,p+1); %minmal value of h
default=struct('intercept',1,'intadjust',0,'alpha',alfa,'h',hdef,'plots',1,'ntrial',dtrial,'classic',0,'Hsets',[]);
list=fieldnames(default);
options=default;
IN=length(list);
i=1;
counter=1;
%Reading optional inputarguments
if nargin > 3
    %
    % placing inputfields in array of strings
    %
    for j=1:nargin-2
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
            error('Both inputarguments alpha and h are provided. Only one is required.')
    end
    %
    % Checking which default parameters have to be changed
    % and keep them in the structure 'options'.
    %
    while counter<=IN 
        index=strmatch(list(counter,:),chklist,'exact');
        if ~isempty(index) % in case of similarity
            for j=1:nargin-3 % searching the index of the accompanying field
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
        if options.h < hmin                                                      
            error(['The LTS must cover at least ' int2str(hmin) ' observations.'])
        elseif options.h > n
            error('h is greater than the number of non-missings and non-infinites.')
        elseif options.h < p
            error(['h should be larger than the dimension ' int2str(p) '.'])
        elseif options.h==0
            options.h=h;
        end
        options.alpha=options.h/n;
    elseif dummy==2
        if options.alpha < 0.5
            options.alpha=0.5;
            mess=sprintf(['Attention : Alpha should be larger than 0.5. \n',...
                    'It is set to 0.5.']);
            disp(mess)
        end
        if options.alpha > 1
            options.alpha=0.75;
            mess=sprintf(['Attention : Alpha should be smaller than 1.\n',...
                    'It is set to 0.75.']);
            disp(mess)
        end
        if options.alpha == 0
            options.alpha=0.75;
        end
        options.h=quanf(options.alpha,n,p);
    end
end
intercept=options.intercept; %if 1 intercept in the model
h=options.h; %number of regular data points on which estimates are based (=[alpha * n])
plots=options.plots; %relevant plots if equal to 1
alfa=options.alpha; %proportion of regular data points
ntrial=options.ntrial; %number of subsets to be taken in the first step
classic=options.classic; %classic least squares regression?
bestobj=inf; %best objective value until 'now' -> goal is to be as small as possible
intadj=options.intadjust; %intercept adjustment needed if set to 1
Hsets = options.Hsets;
if ~isempty(Hsets)
    Hsets_ind = 1;
else
    Hsets_ind = 0;
end

if classic
    res.classic=ols(x,y,'plots',0);
else 
    res.classic=0;
end

if intercept == 1
    dx=dx+[0 1];
    x=cat(2,x,ones(n,1));
end
p=dx(2);
if n < p
   error('Need more observations than variables.');
end

if p > pmax
   error(['The program allows for at most ' int2str(pmax) ' variables.'])
end

rk=rank(x);
if rk < p
   error('x is singular');
end

if h == n 
    res.method='Least Squares Regression.';
    [Q,R]=qr(x,0);
    z=R\(Q'*y);
    raw.coefficients=z;
    residuals=y-x*z;
    raw.res=residuals;
    fitted=x*raw.coefficients;
    raw.fitted=fitted;
    s0=sqrt(sum(residuals.^2)/(n-p));
    if abs(s0) < 1e-7
        weights=abs(residuals)<=1e-7;
        raw.wt=weights;
        raw.scale=0;
        res.scale=0;
        res.coefficients=raw.coefficients;
        raw.objective=0;
    else
        sor=sort(residuals.^2);
        raw.objective=sum(sor(1:h));
        raw.scale=s0;
        weights=abs(residuals/s0)<=norminv(0.9875);
        raw.wt=weights;
        [Q,R]=qr(x(weights==1,:),0);
        z=R\(Q'*y(weights==1));
        res.coefficients=z;
        fitted=x*res.coefficients;
        residuals=y-x*z;
        res.scale=sqrt(sum(weights.*(residuals.^2))/(sum(weights)-1));
        weights=abs(residuals/res.scale)<=norminv(0.9875);
    end
    if intercept
        s1=sum(residuals.^2);
        center=mean(y);
        sh=sum((y-center).^2);
        res.rsquared=1-s1/sh;
    else
        s1=sum(residuals.^2);
        sh=sum(y.^2);
        res.rsquared=1-s1/sh;
    end
    if res.rsquared > 1
        res.rsquared=1;
    elseif res.rsquared < 0
        res.rsquared=0;
    end
    if abs(s0) < 0
        res.method=strvcat(res.method,'An exact fit was found!');
    end
    stdres=residuals/res.scale;
    cutoff.resd=sqrt(chi2inv(0.975,1));
    raw=struct('coefficients',{raw.coefficients},'fitted',{raw.fitted},'res',{raw.res},'scale',{raw.scale},...
        'objective',{raw.objective},'wt',{raw.wt});
    rew=struct('slope',{res.coefficients(1:p)},'int',{0},'fitted',{fitted},'res',{residuals},...
        'scale',{res.scale},'rsquared',{res.rsquared},'h',{h},'alpha',{alfa},'resd', {stdres},...
        'rd',{NaN},'cutoff',{cutoff},'flag',{NaN},'weights',{raw.wt},...
        'method',{res.method},'class',{'LTS'},'classic',{res.classic},'X',{x},'y',{y});
    if intercept
        rew=setfield(rew,'int',res.coefficients(p));
        rew=setfield(rew,'slope',res.coefficients(1:p-1));
    end
    if plots& classic 
        mcdres=mcdcov(x,'h',h,'plots',0);
        if -log(abs(det(mcdres.cov)))/size(data,2)> 50
            res.rd=NaN;
        else 
            res.rd=mcdres.rd;
        end
        cutoff.rd=mcdres.cutoff.rd;  
        cutoff.md=mcdres.cutoff.md;
        flags=abs(stdres)<=cutoff.resd;
        rew=setfield(rew,'rd', res.rd);
        rew=setfield(rew,'flag',flags);
        rew=setfield(rew,'cutoff',cutoff);
        try
            makeplot(rew,'classic',1)
        catch %output must be given even if plots are interrupted 
            %> delete(gcf) to get rid of the menu 
        end
    elseif plots 
        mcdres=mcdcov(x,'h',h,'plots',0);
        if -log(abs(det(mcdres.cov)))/size(x,2)> 50
            res.rd=NaN;
        else 
            res.rd=mcdres.rd;
        end
        cutoff.rd=mcdres.cutoff.rd;  
        cutoff.md=mcdres.cutoff.md;
        flags=abs(stdres)<=cutoff.resd;
        rew=setfield(rew,'rd', res.rd);
        rew=setfield(rew,'cutoff',cutoff);
        rew=setfield(rew,'flag',flags);
        try
            makeplot(rew)
        catch %output must be given even if plots are interrupted 
            %> delete(gcf) to get rid of the menu 
        end
    end
    return
end

if p < 5
    eps=1e-12;
elseif p <= 8
    eps=1e-14;
else
    eps=1e-16;
end

% standardization of the data
xorig=x;
yorig=y;
data=[x y];
if ~intercept
    datamed=repmat(0,1,p+1);
    datamad=median(abs(data)).*1.4826;
    for i=1:p+1
        if abs(datamad(i)) <= eps
            datamad(i)=sum(abs(data(:,i)));
            datamad(i)=(datamad(i)/n)*1.2533;
            if abs(datamad(i)) <= eps;
                error('The MAD of some variable is zero');
            end
        end   
    end
    x=x./repmat(datamad(1:p),n,1);
    y(:,1)=y(:,1)./datamad(p+1);
else
    datamed=median(data);
    datamed(p)=0;
    datamad(p)=1;
    for i=1:p+1
        if i ~= p
            datamad(i)=median(abs(data(:,i)-datamed(i)))*1.4826;
            if abs(datamad(i)) <= eps
                datamad(i)=sum(abs(data(:,i)-datamed(i)));
                datamad(i)=(datamad(i)/n)*1.2533;
                if abs(datamad(i)) <= eps
                    error('The MAD of some variable is zero');
                end
            end
        end
    end
    x=(x-repmat(datamed(1:p),n,1))./repmat(datamad(1:p),n,1);
    y(:,1)=(y(:,1)-datamed(p+1))./datamad(p+1);
end

res.method='Least Trimmed Squares Regression.';
al=0;
teller = zeros(1,n+1);
if Hsets_ind
    csteps = csteps1;
    inplane = NaN;
    fine = 0;
    part =  0;
    final = 1;
    tottimes = 0;
    nsamp = size(Hsets,1);
    obsingroup = n;
else
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
        adjh=floor(group(1)*alfa);
        nsamp=floor(ntrial/ngroup);
        minigr=sum(group);
        obsingroup=fillgroup(n,group,ngroup,seed);
        totgroup=ngroup;
    else
        [part,group,ngroup,adjh,minigr,obsingroup]=deal(0,n,1,h,n,n);
        replow=[50,22,17,15,14,zeros(1,45)];
        if n < replow(p)
            al=1;
            perm=[1:p-1,p-1];                 
            nsamp=nchoosek(n,p);
        else
            al=0;
            nsamp=ntrial;
        end
    end
    
    csteps=csteps1;
    [tottimes,fine,final]=deal(0);
    if part
        bobj1=repmat(inf,ngroup,10); 
        bcoeff1=cell(ngroup,10);
        [bcoeff1{:}]=deal(NaN);
    end
end
bcoeff=cell(1,10);
bobj=repmat(inf,1,10);
[bcoeff{:}]=deal(NaN);
seed=0;
coeffs=repmat(NaN,p,1);
while final ~= 2
    if fine | (~part & final)
        if ~Hsets_ind
            nsamp=10;
        end
        if final
            adjh=h;
            ngroup=1;
            if n*p <= 1e+5
                csteps=csteps3;
            elseif n*p <= 1e+6
                csteps=10-(ceil(n*p/1e+5)-2);
            else
                csteps=1;
            end
            if n > 5000
                nsamp=1;
            end
        else
            adjh=floor(minigr*alfa);
            csteps=csteps2;
        end
    end
    for k=1:ngroup
        for i=1:nsamp
            tottimes=tottimes+1;
            prevobj=0;
            if ~Hsets_ind
                if final
                    if ~isinf(bobj(i))
                        z=bcoeff{i};  	
                    else
                        break
                    end
                elseif fine
                    if ~isinf(bobj1(k,i))
                        z=bcoeff1{k,i};
                    else
                        break
                    end
                else
                    z(1,1)=Inf;
                    while abs(z(1,1)) == Inf
                        if ~part
                            if al
                                k=p;
                                perm(k)=perm(k)+1;
                                while ~(k==1|perm(k) <= (n-(p-k)))
                                    k=k-1;
                                    perm(k)=perm(k)+1;
                                    for j=(k+1):p
                                        perm(j)=perm(j-1)+1;
                                    end	
                                end
                                index=perm;
                            else
                                [index,seed]=randomset(n,p,seed);
                            end
                        else
                            [index,seed]=randomset(group(k),p,seed);
                            index=obsingroup{k}(index);
                        end
                        if p > 1
                            z=x(index,:)\y(index,1);
                            %problems arise whenever the subsample contains
                            %equal x-values. To avoid warnings the tests in line
                            %551 and line 591 are adapted to abs(z(1,1)).
                            %However, in the first run the matrix will still
                            %be singular having coefficients [-inf a b ...] or
                            %[inf inf ...] producing the warning. To avoid
                            %this we turned off the warnings in this
                            %function (line 140).
                        elseif x(index,1) ~= 0
                            z(1,1)=y(index,1)/x(index,1);
                        else
                            z(1,1)=x(index,1);
                        end
                        if al
                            break
                        end
                    end
                end
                if abs(z(1,1)) ~= Inf 
                    if ~part | final
                        residu=y-x*z;
                    elseif ~fine
                        residu=y(obsingroup{k},1)-x(obsingroup{k},:)*z;
                    else
                        residu=y(obsingroup{totgroup+1},1)-x(obsingroup{totgroup+1},:)*z;
                    end        
                    more1=0;
                    more2=0;
                    nmore=200;
                    nmore2=nmore/2;
                    if intadj %intercept adjustment
                        [sortres,sortind]=sort(residu);
                        if ~part %n<600 
                            [center,cover,loc]=mcduni(sortres,obsingroup,adjh,obsingroup-adjh+1,alfa);
                            z(p)=z(p)+center;
                            residu=residu-center;
                        elseif ~fine %n>600, first step
                            [center,cover,loc]=mcduni(sortres,size(obsingroup{k},2),adjh,size(obsingroup{k},2)-adjh+1,alfa);   
                            z(p)=z(p)+center;
                            residu=residu-center;
                        elseif ~final & size(obsingroup{totgroup+1},2)-adjh <= nmore %fine =  merged set
                            [center,cover,loc]=mcduni(sortres,size(obsingroup{totgroup+1},2),adjh,size(obsingroup{totgroup+1},2)-adjh+1,alfa);
                            z(p)=z(p)+center;
                            residu=residu-center;
                        elseif final & n-adjh <= nmore %final = complete data set
                            [center,cover,loc]=mcduni(sortres,n,adjh,n-adjh+1,alfa);
                            z(p)=z(p)+center;
                            residu=residu-center;
                        else   
                            [sortres1,sortind1]=sort(abs(sortres));
                            [sortres2,sortind2]=sort(sortres(sortind1(1:adjh)));
                            further = 1;
                            if final & (sortind1(sortind2(1))+nmore-nmore2+adjh-1 > n  | sortind1(sortind2(1))-nmore2< 1)
                                [center,cover,loc]=mcduni(sortres,n,adjh,n-adjh+1,alfa);
                                z(p)=z(p)+center;
                                residu=residu-center;
                            elseif ~final & fine & (sortind1(sortind2(1))+nmore-nmore2+adjh-1 > size(obsingroup{totgroup+1},2) | sortind1(sortind2(1))-nmore2 < 1)
                                [center,cover,loc]=mcduni(sortres,size(obsingroup{totgroup+1},2),adjh,size(obsingroup{totgroup+1},2)-adjh+1,alfa);
                                z(p)=z(p)+center;
                                residu=residu-center;
                            else   
                                while further
                                    sortres2(1:adjh+nmore)=sortres(sortind1(sortind2(1))-nmore2:sortind1(sortind2(1))+adjh-1+nmore-nmore2);
                                    [center,cover,loc]=mcduni(sortres2,adjh+nmore,adjh,nmore+1,alfa);
                                    if loc == 1 & ~more1
                                        if ~more2
                                            nmore=nmore2;
                                            nmore2=nmore2+nmore2;
                                            more1=1;
                                            if sortind1(sortind2(1))-nmore2 < 1
                                                further=0;
                                            end
                                        else 
                                            further=0;
                                        end
                                    else
                                        if loc == nmore+1 & ~more2
                                            if ~more1
                                                nmore=nmore2;
                                                nmore2=-nmore2;
                                                more2=1;
                                                if final & sortind1(sortind2(1))+nmore-nmore2+adjh-1 > n  
                                                    further=0;
                                                elseif fine & (sortind1(sortind2(1))+nmore-nmore2+adjh-1 > size(obsingroup{totgroup+1},2) | sortind1(sortind2(1))-nmore2<1)
                                                    further=0;
                                                end
                                            else 
                                                further = 0;
                                                
                                            end
                                        else
                                            if loc == 1 & more1
                                                if ~more2
                                                    nmore2=nmore2+100;
                                                    if sortind1(sortind2(1))-nmore2 < 1
                                                        further=0;
                                                    end
                                                else 
                                                    further = 0;
                                                end
                                            else
                                                if loc == nmore+1 & more2
                                                    if ~more1
                                                        nmore2=nmore2+100;
                                                        if final & sortind1(sortind2(1))+nmore-nmore2+adjh-1 > n
                                                            further=0;
                                                        elseif fine & (sortind1(sortind2(1))+nmore-nmore2+adjh-1 > size(obsingroup{totgroup+1},2) | sortind1(sortind2(1))-nmore2<1)
                                                            further=0;
                                                        end
                                                    else 
                                                        further=0;
                                                    end
                                                else
                                                    further=0;
                                                end
                                            end
                                        end
                                    end
                                end
                                z(p)=z(p)+center;
                                residu=residu-center;
                            end
                        end
                    end
                end
            end
            for j=1:csteps %csteps on the subsets
                tottimes=tottimes+1;
                if ~Hsets_ind
                    if z(1,1)~=inf
                        [sortres,sortind]=sort(abs(residu));
                        if fine & ~final
                            sortind=obsingroup{totgroup+1}(sortind);                    
                        elseif part & ~final
                            sortind=obsingroup{k}(sortind);
                        end
                        obs_in_set=sort(sortind(1:adjh));
                        teller(obs_in_set) = teller(obs_in_set) + 1;
                        teller(end) = teller(end) + 1;
                    end
                else
                    obs_in_set = Hsets(i,:);
                end
                
                if Hsets_ind|(~Hsets_ind & z(1,1)~=inf)
                    [Q,R]=qr(x(obs_in_set,:),0);
                    z=R\(Q'*y(obs_in_set,1));
                    if ~part | final
                        residu=y-x*z;   
                    elseif ~fine
                        residu=y(obsingroup{k},1)-x(obsingroup{k},:)*z;
                    else
                        residu=y(obsingroup{totgroup+1},1)-x(obsingroup{totgroup+1},:)*z;
                    end
                    more1=0;
                    more2=0;
                    nmore=200;
                    nmore2=nmore/2;
                    if intadj %intercept adjustment
                        [sortres,sortind]=sort(residu);
                        if ~part
                            [center,cover,loc]=mcduni(sortres,obsingroup,adjh,obsingroup-adjh+1,alfa);
                            z(p)=z(p)+center;
                            residu=residu-center;
                        elseif ~fine
                            [center,cover,loc]=mcduni(sortres,size(obsingroup{k},2),adjh,size(obsingroup{k},2)-adjh+1,alfa);   
                            z(p)=z(p)+center;
                            residu=residu-center;
                        elseif ~final & size(obsingroup{totgroup+1},2)-adjh <= nmore
                            [center,cover,loc]=mcduni(sortres,size(obsingroup{totgroup+1},2),adjh,size(obsingroup{totgroup+1},2)-adjh+1,alfa);
                            z(p)=z(p)+center;
                            residu=residu-center;
                        elseif final & n-adjh <= nmore
                            [center,cover,loc]=mcduni(sortres,n,adjh,n-adjh+1,alfa);
                            z(p)=z(p)+center;
                            residu=residu-center;
                        else   
                            [sortres1,sortind1]=sort(abs(sortres));
                            [sortres2,sortind2]=sort(sortres(sortind1(1:adjh)));
                            further = 1;
                            if final & (sortind1(sortind2(1))+nmore-nmore2+adjh-1 > n  | sortind1(sortind2(1))-nmore2< 1)
                                [center,cover,loc]=mcduni(sortres,n,adjh,n-adjh+1,alfa);
                                z(p)=z(p)+center;
                                residu=residu-center;
                            elseif ~final & fine & (sortind1(sortind2(1))+nmore-nmore2+adjh-1 > size(obsingroup{totgroup+1},2) | sortind1(sortind2(1))-nmore2 < 1)
                                [center,cover,loc]=mcduni(sortres,size(obsingroup{totgroup+1},2),adjh,size(obsingroup{totgroup+1},2)-adjh+1,alfa);
                                z(p)=z(p)+center;
                                residu=residu-center;
                            else   
                                while further
                                    sortres2(1:adjh+nmore)=sortres(sortind1(sortind2(1))-nmore2:sortind1(sortind2(1))+adjh-1+nmore-nmore2);
                                    [center,cover,loc]=mcduni(sortres2,adjh+nmore,adjh,nmore+1,alfa);
                                    if loc == 1 & ~more1
                                        if ~more2
                                            nmore=nmore2;
                                            nmore2=nmore2+nmore2;
                                            more1=1;
                                            if sortind1(sortind2(1))-nmore2 < 1
                                                further=0;
                                            end
                                        else 
                                            further=0;
                                        end
                                    else
                                        if loc == nmore+1 & ~more2
                                            if ~more1
                                                nmore=nmore2;
                                                nmore2=-nmore2;
                                                more2=1;
                                                if final & sortind1(sortind2(1))+nmore-nmore2+adjh-1 > n  
                                                    further=0;
                                                elseif fine & (sortind1(sortind2(1))+nmore-nmore2+adjh-1 > size(obsingroup{totgroup+1},2) | sortind1(sortind2(1))-nmore2<1)
                                                    further=0;
                                                end
                                            else further=0; 
                                            end
                                        else
                                            if loc == 1 & more1
                                                if ~more2
                                                    nmore2=nmore2+100;
                                                    if sortind1(sortind2(1))-nmore2 < 1
                                                        further=0;
                                                    end
                                                else further = 0; 
                                                end
                                            else
                                                if loc == nmore+1 & more2
                                                    if ~more1
                                                        nmore2=nmore2+100;
                                                        if final & sortind1(sortind2(1))+nmore-nmore2+adjh-1 > n
                                                            further=0;
                                                        elseif fine & (sortind1(sortind2(1))+nmore-nmore2+adjh-1 > size(obsingroup{totgroup+1},2) | sortind1(sortind2(1))-nmore2<1)
                                                            further=0;
                                                        end
                                                    else further=0; 
                                                    end
                                                else
                                                    further=0;
                                                end
                                            end
                                        end
                                    end
                                end
                                z(p)=z(p)+center;
                                residu=residu-center;
                            end
                        end
                    end
                    sor=sort(abs(residu));
                    obj=sum(sor(1:adjh).^2); %objective function value after the iteration
                    if j >= 2 & obj == prevobj
                        break;
                    end
                    prevobj=obj;
                end %end final
            end
            if ~final 
                if fine |~part %merged set or n<600
                    if obj < max(bobj)		
                        [bcoeff,bobj]=insertion(bcoeff,bobj,z,obj,1,eps);
                    end
                else
                    if obj < max(bobj1(k,:))	
                        [bcoeff1,bobj1]=insertion(bcoeff1,bobj1,z,obj,k,eps);
                    end
                end
            end
            if final & obj < bestobj
                bestset=obs_in_set;
                bestobj=obj;
                coeffs=z;
            end 
        end    %end for nsamp
    end %end for ngroups
    if part & ~fine
        fine = 1;
    elseif (part & fine & ~final) | (~part & ~final)
        final = 1;
    else
        final = 2;
    end
end %end while,  so final = 2

if p <= 1
    coeffs(1)=coeffs(1)*datamad(p+1)/datamad(1);
else
    coeffs(1:p-1)=(coeffs(1:p-1)*datamad(p+1))'./datamad(1:p-1);
    if ~intercept
        coeffs(p)=coeffs(p)*datamad(p+1)/datamad(p);
    else
        coeffs(p)=coeffs(p)*datamad(p+1);
        coeffs(p)=coeffs(p)-sum(coeffs(1:p-1)'.*datamed(1:p-1));
        coeffs(p)=coeffs(p)+datamed(p+1);
    end      
end
bestobj=bestobj*(datamad(p+1)^2);
x=xorig;
y=yorig;
raw.coefficients=coeffs;
raw.objective=bestobj;
fitted=x*coeffs;
raw.fitted=fitted;
residuals=y-fitted;
raw.residuals=residuals;
sor=sort(residuals.^2);
factor=rawconsfactorlts(h,n);
sh0=sqrt((1/h)*sum(sor(1:h)));
s0=sh0*factor;
cutoff.resd=sqrt(chi2inv(0.975,1));
if abs(s0) < 1e-7
   weights=abs(residuals)<=1e-7;
   raw.wt=weights;
   raw.scale=0;
   res.scale=0;
   res.coefficients=raw.coefficients;
   raw.objective=0;
else
    raw.scale=s0;
    m=2*n/asvarscalekwad(h,n);
    quantile=tinv(0.9875,m);
    weights=abs(residuals/s0)<=quantile;
    raw.wt=weights;
    [Q,R]=qr(x(weights==1,:),0);
    z=R\(Q'*y(weights==1));
    res.coefficients=z;
    fitted=x*res.coefficients;
    residuals=y-fitted;
    res.scale=sqrt(sum(weights.*residuals.^2)/(sum(weights)-1));
    s0=res.scale;
    weights=abs(residuals/res.scale)<=cutoff.resd;
end
res.flag=weights;

res.Hsubsets.Hopt = bestset;
[telobs,indobs] = greatsort(teller(1:(end - 1)));
res.Hsubsets.Hfreq = indobs(1:(h));
if size(res.Hsubsets.Hfreq,2) == 1
    res.Hsubsets.Hfreq = res.Hsubsets.Hfreq';
end

if intercept
    yw=y(raw.wt==1);
    cyw=mcenter(yw);
    sres=sum(residuals(raw.wt==1).^2);
    cwy2=sum(cyw.^2);
    res.rsquared=1-sres/cwy2;
else
    sor=sort(residuals.^2);
    s1=sum(sor(1:h));
    sor=sort(y.^2);
    sh=sum(sor(1:h));
    res.rsquared=1-(s1/sh);
end
if res.rsquared > 1
   res.rsquared=1;
elseif res.rsquared < 0
    res.rsquared=0;
end

if abs(s0) < 1e-7
    res.method=strvcat(res.method,'An exact fit was found!');
    res.Hsubsets.Hopt=1:n;
    res.Hsubsets.Hfreq=1:n;
    disp('Exact fit was encountered')
end   

if ~intercept
    data=x;
else 
    data=x(:,1:p-1);
end
%calculating residual distances : in case of a univariate analysis they are
%equal to the standardized residuals.
stdres=residuals/res.scale;
cutoff.resd=sqrt(chi2inv(0.975,1));

%assigning ouput
raw=struct('coefficients',{raw.coefficients},'fitted',{raw.fitted},'res',{raw.residuals},'scale',{raw.scale},...
    'objective',{raw.objective},'wt',{raw.wt});
rew=struct('slope',{res.coefficients(1:p)},'int',{0},'fitted',{fitted},'res',{residuals},...
    'scale',{res.scale},'rsquared',{res.rsquared},'h',{h},'Hsubsets',{res.Hsubsets},'alpha',{alfa},'rd',{0},'resd', {stdres},...
    'cutoff',{cutoff},'flag',{res.flag},...
    'method',{res.method},'class',{'LTS'},'classic',{res.classic},'X',{data},'y',{y});

if intercept
    rew=setfield(rew,'int',res.coefficients(p));
    rew=setfield(rew,'slope',res.coefficients(1:p-1));
end

if isfield(rew,'X') & ((size(x,2)-intercept)~=1 | size(y,2)~=1)
    rew=rmfield(rew,{'X','y'});
end
warning on
if plots & classic
    mcdres=mcdcov(data,'h',h,'plots',0);
    if -log(abs(det(mcdres.cov)))/size(data,2)> 50
        res.rd=NaN;
    else
        res.rd=mcdres.rd;
    end
    cutoff.rd=mcdres.cutoff.rd; 
    cutoff.md=mcdres.cutoff.md;
    rew=setfield(rew,'cutoff',cutoff);
    rew=setfield(rew,'rd', res.rd);
    try
        makeplot(rew,'classic',1)
    catch %output must be given even if plots are interrupted 
        %> delete(gcf) to get rid of the menu 
    end
elseif plots
    mcdres=mcdcov(data,'h',h,'plots',0);
    if -log(abs(det(mcdres.cov)))/size(data,2)> 50
        res.rd=NaN;
    else
        res.rd=mcdres.rd;
    end
    cutoff.rd=mcdres.cutoff.rd; 
    cutoff.md=mcdres.cutoff.md;
    rew=setfield(rew,'cutoff',cutoff);
    rew=setfield(rew,'rd', res.rd);
    try
        makeplot(rew)
    catch %output must be given even if plots are interrupted 
        %> delete(gcf) to get rid of the menu 
    end
end
if rew.rd==0
    rew=rmfield(rew,'rd');
end
% --------------------------------------------------------------------

function obsingroup = fillgroup(n,group,ngroup,seed)

% Creates the subdatasets.

obsingroup=cell(1,ngroup+1);
jndex=0;
for k = 1:ngroup
   for m = 1:group(k)
      [random,seed]=uniran(seed);
      ran=floor(random*(n-jndex)+1);
      jndex=jndex+1;
      if jndex == 1
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

% --------------------------------------------------------------------

function [ranset,seed] = randomset(tot,nel,seed)

for j = 1:nel
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


% --------------------------------------------------------------------

function [output] = replow(k,pmax)

replow=[500,50,22,17,15,14];
help=zeros(1,pmax-5);
replow=[replow help];
output=replow(k);

% --------------------------------------------------------------------

function [initmean,initcov,iloc]=mcduni(y,ncas,h,len,alfa)

% The exact MCD algorithm for the univariate case. 

y=sort(y);
ay(1)=sum(y(1:h));
factor=rawconsfactormcd(h,ncas);

for samp=2:len
   ay(samp)=ay(samp-1)-y(samp-1)+y(samp+h-1);
end
ay2=ay.^2/h;
sq(1)=sum(y(1:h).^2)-ay2(1);
for samp=2:len
   sq(samp)=sq(samp-1)-y(samp-1)^2+y(samp+h-1)^2-ay2(samp)+ay2(samp-1);
end
sqmin=min(sq);
ii=find(sq==sqmin);
ndup=length(ii);
slutn(1:ndup)=ay(ii);
initmean=slutn(floor((ndup+1)/2))/h;
initcov=factor^2*sqmin/h;
iloc=ii(1);

% -----------------------------------------------------------------------------

function [bestmean,bobj]=insertion(bestmean,bobj,z,obj,row,eps)

insert=1;
equ=find(obj==bobj(row,:));
for j=equ
   if (z==bestmean{row,j}) 
      insert=0;
   end
end
if insert 
   ins=min(find(obj < bobj(row,:))); 
   if ins==10
      bestmean{row,ins}=z;
      bobj(row,ins)=obj;
   else
      [bestmean{row,ins+1:10}]=deal(bestmean{row,ins:9});
      bestmean{row,ins}=z;
      bobj(row,ins+1:10)=bobj(row,ins:9);
      bobj(row,ins)=obj;
   end
end

% --------------------------------------------------------------------------------

function quan=quanf(alfa,n,rk)

quan=floor(2*floor((n+rk+1)/2)-n+2*(n-floor((n+rk+1)/2))*alfa);


%---------------------------------------------------------------

function rawconsfacmcd=rawconsfactormcd(quan,n)

qalpha=chi2inv(quan/n,1);
calphainvers=gamcdf(qalpha/2,1/2+1)/(quan/n);
calpha=1/calphainvers;
rawconsfacmcd=calpha;


%-------------------------------------------------------------

function rawconsfaclts=rawconsfactorlts(quan,n)

rawconsfaclts=(1/sqrt(1-((2*n)/(quan*(1/norminv((quan+n)/(2*n)))))*...
		normpdf(1/(1/(norminv((quan+n)/(2*n)))))));

%--------------------------------------------------------------

function asvar=asvarscalekwad(quan,n)

alfa=quan/n;
alfa=1-alfa;
qalfa=chi2inv(1-alfa,1);
c2=gamcdf(qalfa/2,1/2+1);
c1=1/c2;
c3=3*gamcdf(qalfa/2,1/2+2);
asvar=qalfa*(1-alfa)-c2;
asvar=asvar^2;
asvar=(c3-2*qalfa*c2+(1-alfa)*(qalfa^2))-asvar;
asvar=c1^2*asvar;
%--------------------------------------------------------------




