function result=rpcr(x,y,varargin)

%RPCR is a 'Robust Principal Components Regression' method based on ROBPCA.
% It can be applied to both low and high-dimensional predictor variables x,
% and to one or multiple response variables y. It is resistant to outliers 
% in the data. First, a robust principal components method is applied to the 
% predictor variables x (see robpca.m). Then a robust regression is performed. 
% For univariate y, the LTS regression is used (see ltsregres.m). When there are 
% several response variables, the MCD-regression method is applied (see mcdregres.m).
%
% The RPCR method is described in: 
%    Hubert, M., Verboven, S. (2003),
%    "A robust PCR method for high-dimensional regressors",
%    Journal of Chemometrics, 17, 438-452.
%
% To select the number of components in the regression model, a robust RMSECV (root mean squared
% error of cross validation) curve is drawn, based on a fast algorithm for
% cross-validation. This approach is described in:
%
%    Engelen, S., Hubert, M. (2005),
%    "Fast model selection for robust calibration methods",
%    Analytica Chimica Acta, 544, 219-228.
%
% Required input arguments:
%            x : Data matrix of the explanatory variables
%                (n observations in rows, p variables in columns)
%            y : Data matrix of the response variables
%                (n observations in rows, q variables in columns)
%
% Optional input arguments:
%            k : Number of principal components to be used in the PCA step 
%                (default = min(rank(x),kmax)). If k is not specified,
%                it can be selected using the option 'rmsecv'. 
%         kmax : Maximal number of principal components to be used (default = 10).
%                If k is provided, kmax does not need to be specified, unless k is larger
%                than 10.                 
%        alpha : (1-alpha) measures the fraction of outliers the algorithm should 
%                resist. Any value between 0.5 and 1 may be specified. (default = 0.75)
%            h : (n-h+1) measures the number of observations the algorithm should 
%                resist. Any value between n/2 and n may be specified. (default = 0.75*n)
%                Alpha and h may not both be specified.
%       rmsecv : If equal to zero and k is not specified, a robust R-squared curve
%                is plotted and an optimal k value can be chosen. (default)
%                If equal to one and k is not specified, a robust component selection-curve is plotted and
%                an optimal k value can be chosen. This curve computes a combination of
%                the robust cross-validated mean squared error and the robust residual
%                sum of squares for k=1 to kmax. Rmsecv and k may not both
%                be specified.
%        rmsep : If equal to one, the robust RMSEP-value (root mean squared error of
%                prediction) for the model with k components.                
%                (default = 0). This value is automatically given if rmsecv = 1.
%    intadjust : If equal to one, the intercept adjustment for the
%                LTS-regression will be calculated (default = 0). See ltsregres.m for
%                details.
%        plots : If equal to one, a menu is shown which allows to draw several plots,
%                such as a score outlier map and a regression outlier map. (default)
%                If the input argument 'classic' is equal to one, the classical
%                plots are drawn as well.
%                If 'plots' is equal to zero, all plots are suppressed.
%                See also makeplot.m
%        labsd : The 'labsd' observations with largest score distance are
%                labeled on the diagnostic plots. (default = 3)
%        labod : The 'labod' observations with largest orthogonal distance are
%                labeled on the score outlier map. (default = 3)    
%      labresd : The 'labresd' observations with largest residual distance are
%                labeled on the regression outlier map. (default = 3)   
%      classic : If equal to one, the classical PCR analysis will be performed as well
%                (see also cpcr.m). (default = 0)
% 
%
% I/O: result=rpcr(x,y,'k',0,'kmax',10,'alpha',0.75,'h',h,'rmsecv',0,'rmsep',0,...
%             'plots',1,'labsd',3,'labod',3,'labresd',3,'classic',0);
%  The user should only give the input arguments that have to change their default value.
%  The name of the input arguments needs to be followed by their value.
%  The order of the input arguments is of no importance.
%   
% Example: result=rpcr(x,y,'alpha',0.65,'k',5,'plots',0,'classic',1);
%          result=rpcr(x,y,'kmax',5,'labresd',7);
%
% The output of RPCR is a structure containing:
%
%   result.slope     : Robust slope
%   result.int       : Robust intercept 
%   result.fitted    : Robust prediction vector
%   result.res       : Robust residuals
%   result.cov       : Estimated variance-covariance matrix of the residuals  
%   result.rsquared  : Robust R-squared value for the optimal k. 
%   result.rcs       : Robust Component Selection Criterion:
%                      This is a matrix with kmax columns. The first row contains the approximate
%                      R-squared values (for k=1,...,kmax) and the second row the square root of the weighted
%                      residuals sum of squares. This is equal to the RCS-value with 
%                      gamma = 0 (see Engelen and Hubert (2005) for the definition). 
%                      If the input argument rmsecv = 1, the third row contains the RCS-value sfor
%                      gamma = 0.5 and the fourth for gamma = 1. The last one is equal to the robust 
%                      cross-validated RMSE values. 
%                      Note that all the entries in this matrix depend on the choice of kmax. 
%   result.rmsep     : Robust RMSEP value 
%   result.k         : Number of principal components
%   result.h         : The quantile h used throughout the algorithm
%   result.sd        : Robust score distances within the robust PCA subspace
%   result.od        : Robust orthogonal distances to the robust PCA subspace 
%   result.resd      : Residual distances (when there are several response variables).
%                      If univariate regression is performed, it contains the standardized residuals.
%   result.cutoff    : Cutoff values for the score, orthogonal and residual distances
%   result.flag      : The observations whose score distance is larger than 
%                      'result.cutoff.sd' receive a flag 'result.flag.sd' equal
%                      to zero (good leverage points). Otherwise 'result.flag.sd'
%                      is equal to one. 
%                      The components 'result.flag.od' and 'result.flag.resd' are
%                      defined analogously, and determine the orthogonal outliers, 
%                      resp. the bad leverage points/vertical outliers. 
%                      The observations with 'result.flag.od' and 'result.flag.resd'
%                      equal to zero, can be considered as calibration outliers and receive
%                      'result.flag.all' equal to zero. The regular observations and the good leverage
%                      points have 'result.flag.all' equal to one.
%   result.class     : 'RPCR'
%   result.classic   : If the input argument 'classic' is equal to one, this structure
%                      contains results of the classical PCR analysis (see also cpcr.m). 
%   result.robpca    : Full output of the robust PCA analysis (see robpca.m) 
%   result.lts       : If there is one response variable: full output of the LTS regression.
%                      (see ltsregres.m)
%   result.mcdreg    : If there are several response variables: full output of the MCD regression.
%                      (see mcdregres.m)
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Sabine Verboven
% Created on: 31/10/2000
% Last Update: 09/04/2004,  03/07/2006
% Last Revision: 04/08/2006

%
% checking input
%
if rem(nargin,2)~=0
    error('Number of inputarguments must be even!');
end
%
% initialization with defaults
%
counter=1;
[n,p]=size(x);
r=rank(x);
[n2,q]=size(y);
if n~=n2
    error('The response variables and the predictor variables have a different number of observations.')
end
kmax=min([10,floor(n/2),r]);
alfa=0.75;
k=0;
h=floor(2*floor((n+kmax+2)/2)-n+2*(n-floor((n+kmax+2)/2))*alfa);
default=struct('intadjust', 0, 'alpha',alfa,'h',h,'k',k,'kmax',kmax,'plots',1,...
    'rmsecv',0,'classic',0,'labsd',3,'labod',3,'labresd',3,'rmsep',0);
list=fieldnames(default);
options=default;
IN=length(list);
i=1;
%
if nargin>2
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
            options.alpha=alfa;
            options.h=h;
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
    options.h=floor(options.h);
    options.kmax=floor(options.kmax);
    options.k=floor(options.k);
    options.labsd=max(0,min(floor(options.labsd),n));
    options.labod=max(0,min(floor(options.labod),n));
    options.labresd=max(0,min(floor(options.labresd),n));
    kmax=min([options.kmax,floor(n/2),r]);
    dummyh = strcmp(chklist,'h');
    dummykmax = strcmp(chklist,'kmax');
    if all(dummyh == 0) && any(dummykmax)
        h = floor(2*floor((n+kmax+1)/2)-n+2*(n-floor((n+kmax+1)/2))*alfa);
    end
    if options.k < 0 
        options.k=0;
    elseif options.k > kmax 
        mess=sprintf(['Attention (rpcr.m): The required number of principal components, k = ',num2str(options.k)...
            ,'\n is larger than kmax= ',num2str(kmax),'; k is set to ',num2str(kmax)]);
        disp(mess)
        options.k=kmax;
    end
    if dummy==1% checking inputvariable h
        if options.h-floor(options.h)~=0
            mess=sprintf(['Attention (rpcr.m): h must be an integer. \n']);
            disp(mess)
        end
        if options.k==0
            if options.h<floor((n+kmax+1)/2 )
                options.h = quanf(0.5,n,kmax);
                mess=sprintf(['Attention (rpcr.m): h should be larger than (n+2)/2.\n',...
                        'It is set to its minimum value ',num2str(options.h)]);
                disp(mess)
            end
            options.alpha=options.h/n;
        else
            if options.h<floor((n+options.k+1)/2) % checking inputvariable h
                options.h = quanf(0.5,n,k);
                mess=sprintf(['Attention (rpcr.m): h should be larger than (n+2)/2.\n',...
                        'It is set to its minimum value ',num2str(options.h)]);
                disp(mess)
            end
            options.alpha=options.h/n;
        end
        if options.h>n
            options.alpha=0.75;
            if options.k==0
                options.h=floor(2*floor((n+kmax+2)/2)-n+2*(n-floor((n+kmax+2)/2))*options.alpha);
            else
                options.h=floor(2*floor((n+options.k+1)/2)-n+2*(n-floor((n+options.k+1)/2))*options.alpha);%k+1 ipv k?
            end
            mess=sprintf(['Attention (rpcr.m): h should be smaller than n. \n',...
                    'It is set to its default value ',num2str(options.h)]);
            disp(mess)
        end
    elseif dummy==2
        if options.alpha < 0.5
            options.alpha=0.5;
            mess=sprintf(['Attention (rpcr.m): Alpha should be larger than 0.5. \n',...
                    'It is set to 0.5.']);
            disp(mess)
        end
        if options.alpha > 1
            options.alpha=0.75;
            mess=sprintf(['Attention (rpcr.m): Alpha should be smaller than 1.\n',...
                    'It is set to 0.75.']);
            disp(mess)
        end
        if options.k==0
            options.h=floor(2*floor((n+kmax+2)/2)-n+2*(n-floor((n+kmax+2)/2))*options.alpha);
        else
            options.h=floor(2*floor((n+options.k+1)/2)-n+2*(n-floor((n+options.k+1)/2))*options.alpha);%k+1 ipv k?
        end
    end
end
%
%
% MAIN PART
%
% If y is univariate, perform LTS regression, else MCD regression
if q==1 
    options.method='lts';
else
    options.method='mcd';
end
% Default value for the number of latent variables k
if n<kmax+2
     error('You need more data points to compute this estimator.')
end
%  
% initializing and checking k
if options.k==0 % no optimal number of latent variables given by user
    kfixed = 0;
    if options.rmsecv==0 % no cross validation asked -> calculate R2
        [R2,final]=rsquared(x,y,kmax,'RPCR',options.h);
        rss = final.rss;
        options.k=final.k;
        outrobpca=robpca(x,'h',options.h,'k',options.k,'kmax',kmax,'plots',0,'classic',options.classic);
        options.k=outrobpca.k;
        final=regression(x,y,options.k,outrobpca,options);
        if options.rmsep==1 % the rmsep is asked, without prior knowledge
            rmse=rrmse(x,y,options.h,kmax,'RPCR',0,options.k);
            final.rmsep=rmse.rmsep;
        end
        final.robpca=outrobpca;
        final.R2 = R2;
        final.rss = rss;
    else % RMSE with cross validation
        crosscv=rrmse(x,y,options.h,kmax,'RPCR'); 
        options.k=crosscv.k;
        rmsecv=crosscv.rmsecv;
        pred=rrmse(x,y,options.h,kmax,'RPCR',0,options.k,crosscv.weight,crosscv.res); %calculate rmsep with rmsecv weights and residuals
        outrobpca=robpca(x,'h',options.h,'k',options.k,'kmax',kmax,'plots',0,'classic',options.classic);
        options.k=outrobpca.k;
        final=regression(x,y,options.k,outrobpca,options);
        final.robpca=outrobpca;
        final.rmsecv=rmsecv;
        final.rmsep=pred.rmsep;
        final.R2 = crosscv.R2;
        final.rss = crosscv.rss;
    end
else % optimal number of latent variables is given by user
    if options.rmsecv
        error(['Both RMSECV and k were given.', ...
            'Please rerun your analysis with one of these inputs.  (see help file)'])
    end
    kfixed = 1;
    k=min(r,min(options.k,kmax));
    options.k=k;
    outrobpca=robpca(x,'h',options.h,'k',options.k,'kmax',kmax,'plots',0,'classic',options.classic);
    options.k=outrobpca.k;
    final=regression(x,y,options.k,outrobpca,options);
    if options.rmsep==1 % the rmsep is asked, without prior knowledge
        rmse=rrmse(x,y,options.h,kmax,'RPCR',0,options.k);
        final.rmsep=rmse.rmsep;
    end
    final.robpca=outrobpca;
end

final.h=options.h;
final.k=options.k;
final.alpha=options.alpha;

% scores-distances
final.sd=final.robpca.sd;
quan=chi2inv(0.975,final.k);
final.cutoff.sd=quan^0.5;

% orthogonal distances 
final.od=final.robpca.od;
final.cutoff.od=final.robpca.cutoff.od;
final.x=x;
final.y=y;

if options.classic==1
    final.classic=cpcr(x,y,'k',final.k,'plots',0);
else
    final.classic=0;
end

if strmatch(options.method,'mcd','exact') 
    final.resd=final.mcdreg.resd;
    final.cutoff.resd=final.mcdreg.cutoff.resd;
else
    final.resd=final.lts.res/final.lts.scale;
    final.cutoff.resd=sqrt(chi2inv(0.975,1));
end
final.class='RPCR';

%Computing flags
final.flag.od=final.od<=final.cutoff.od;
final.flag.resd=abs(final.resd)<=final.cutoff.resd;
final.flag.all=(final.flag.od & final.flag.resd);

% Assigning output
if kfixed == 0
    final.rcs = [final.R2;sqrt(final.rss)];
else
    final.rcs = 0;
end

if options.rmsecv~=0
    gammahalf = 0.5*sqrt(final.rss) + 0.5*final.rmsecv;
    final.rcs = [final.rcs;gammahalf;final.rmsecv];
    rmsecv_ind = 1;
else
    rmsecv_ind = 0;
end
if options.rmsep~=1 && rmsecv_ind == 0
    final.rmsep=0;
end

if isfield(final,'lts')
    result=struct('slope',{final.slope}, 'int',{final.int}, 'fitted',{final.fitted},'res',{final.res},...
        'cov',{final.cov},'rsquared',{final.rsquared},'rcs',{final.rcs},'rmsep',...
        {final.rmsep},'k',{final.k},'alpha',{final.alpha},'h',{final.h},'sd', {final.sd},'od',{final.od},'resd',{final.resd},...
        'cutoff',{final.cutoff},'flag',{final.flag},'class',{final.class},'classic',{final.classic},...
        'robpca',{final.robpca},'lts',{final.lts});
else
    result=struct('slope',{final.slope}, 'int',{final.int}, 'fitted',{final.fitted},'res',{final.res},...
        'cov',{final.cov},'rsquared',{final.rsquared},'rcs',{final.rcs},'rmsep',{final.rmsep},...
        'k',{final.k},'alpha',{final.alpha},'h',{final.h},'sd', {final.sd},'od',{final.od},'resd',{final.resd},...
        'cutoff',{final.cutoff},'flag',{final.flag},'class',{final.class},'classic',{final.classic},...
        'robpca',{final.robpca},'mcdreg',{final.mcdreg});
end
if result.rcs == 0
    result = rmfield(result,'rcs');
end
if result.rmsep==0 
    result=rmfield(result,'rmsep');
end

% Plots
try
    if options.plots && options.classic
        makeplot(result,'classic',1,'labsd',options.labsd,'labod',options.labod,'labresd',options.labresd)
    elseif options.plots
        makeplot(result,'labsd',options.labsd,'labod',options.labod,'labresd',options.labresd) 
    end
catch %output must be given even if plots are interrupted 
    %> delete(gcf) to get rid of the menu 
end
%---------------------------------
function out=regression(x,y,d,pre_out,options)

switch options.method
case 'mcd' % MCD regression for multivariate responses
    out.mcdreg=mcdregres(pre_out.T(:,1:d),y,'k',options.k,'h',options.h,'plots',0);
    if ~isstruct(out.mcdreg)
        return
    end
    %%coefficients in the original space
    out.slope=pre_out.P(:,1:d)*out.mcdreg.slope;
    out.int=out.mcdreg.int-pre_out.M*pre_out.P(:,1:d)*out.mcdreg.slope;
    out.fitted=[pre_out.T(:,1:d) ones(size(pre_out.T(:,1:d),1),1)]*[out.mcdreg.slope; out.mcdreg.int]; 
    out.res=y-out.fitted;
    out.cov=out.mcdreg.cov; % covariance matrix of residuals (capital sigma)
    out.name='Multivariate robust MCD-regression';  
    out.rsquared=out.mcdreg.rsquared;
case 'lts' 
    [out.lts,raw]=ltsregres(pre_out.T(:,1:d),y,'plots',0,'h',options.h,'intadjust',options.intadjust);
    %%coefficients in the original space;
    out.lts.weights=raw.wt;
    out.slope=pre_out.P(:,1:d)*out.lts.slope;
    out.int=out.lts.int-pre_out.M*pre_out.P(:,1:d)*out.lts.slope;
    out.fitted=[pre_out.T(:,1:d) ones(size(pre_out.T(:,1:d),1),1)]*[out.lts.slope; out.lts.int];
    out.res=y-out.fitted; 
    out.cov=out.lts.scale.^2;
    if out.cov==0
        mess=sprintf(['Attention (rpcr.m): No standardized residuals could be calculated \n',...
            'because their robust scale is zero.']);
        disp(mess)
        out.stdres=NaN;
    else
        out.stdres=out.res/out.lts.scale;
    end
    out.rsquared=out.lts.rsquared;
    out.name='Robust LTS regression';
end

%-----------------------------------------------------------------------------------------
function quan=quanf(alfa,n,rk)

quan=floor(2*floor((n+rk+1)/2)-n+2*(n-floor((n+rk+1)/2))*alfa);





