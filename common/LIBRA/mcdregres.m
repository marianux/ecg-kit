function result=mcdregres(x,y,varargin)

%MCDREGRES is a robust multivariate regression method. It can handle multiple
% response variables. The estimates are based on the robust MCD estimator of 
% location and scatter (see mcdcov.m). The explanatory variables should be 
% low-dimensional, otherwise robust principal component regression (rpcr.m) 
% or robust partial least squares (rsimpls.m) should be applied.
%
% The MCD regression method is described in 
%    Rousseeuw, P.J., Van Aelst, S., Van Driessen, K, Agullo, J. (2004),
%    "Robust multivariate regression", Technometrics, 46, pp 293-305.
%
% Required input arguments:
%            x : Data matrix of the explanatory variables
%                (n observations in rows, p variables in columns)
%            y : Data matrix of the response variables
%                (n observations in rows, q variables in columns)
%
% Optional input arguments: 
%        alpha : (1-alpha) measures the amount of contamination the algorithm should 
%                resist. Any value between 0.5 and 1 may be specified. (default = 0.75)
%            h : The quantile of observations whose covariance determinant will 
%                be minimized.  Any value between n/2 and n may be specified.
%                The default value is 0.75*n.
%       ntrial : The number of random trial subsamples that are drawn for 
%                large datasets. (default = 500)
%        plots : If equal to one, a menu is shown which allows to draw a regression
%                outlier map. (default)
%                If the input argument 'classic' is equal to one, the classical
%                plot is drawn as well.
%                If 'plots' is equal to zero, all plots are suppressed.
%                See also makeplot.m
%      classic : If equal to one, classical multivariate linear regression 
%                is performed as well, see mlr.m. (default = 0)
%
% Input arguments for advanced users:
%     Hsets : Instead of random trial h-subsets (default, Hsets = []), Hsets makes it possible to give certain
%             h-subsets as input. Hsets is a matrix that contains the indices of the observations of one
%             h-subset as a row.
%
% I/O: result=mcdregres(x,y,'alpha',0.75,'ntrial',500,'plots',1,'classic',0);
%  The user should only give the input arguments that have to change their default value.
%  The name of the input arguments needs to be followed by their value.
%  The order of the input arguments is of no importance.
%   
% Example: result=mcdregres(x,y,'plots',0,'alpha',0.70)
%
% The output is a structure which contains
%   result.slope     : Robust slope (matrix)
%   result.int       : Robust intercept (vector)
%   result.fitted    : Robust prediction matrix
%   result.res       : Robust residuals
%   result.cov       : Estimated variance-covariance matrix of the residuals 
%   result.rsquared  : Robust R-squared value
%   result.h         : The quantile h used throughout the algorithm
%   result.Hsubsets  : A structure that contains Hopt and Hfreq:
%                        Hopt  : The subset of h points whose covariance matrix has minimal determinant, 
%                                ordered following increasing robust distances.
%                        Hfreq : The subset of h points which are the most frequently selected during the whole
%                                algorithm.
%   result.rd        : Robust scores distances in x-space
%   result.resd      : Residual distances (when there are several response variables).
%                      If univariate regression is performed, it contains the standardized residuals.
%   result.cutoff    : Cutoff values for the score and residual distances
%   result.weights   : The observations with weight one are used in the reweighting, 
%                      the other observations have zero weight.
%   result.flag      : The observations whose residual distance is larger than result.cutoff.resd
%                      (bad leverage points/vertical outliers) can be considered as outliers and receive 
%                      a flag equal to zero.
%                      The regular observations, including the good leverage points, 
%                      receive a flag 1.
%   result.class     : 'MCDREG'
%   result.classic   : If the input argument 'classic' is equal to one, this structure
%                      contains results of classical multivariate regression (see also mlr.m). 
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Original S-PLUS code by Katrien Vandriessen, implemented in MATLAB by Sabine Verboven 
% Version date : 12/02/2004
% Last update: 04/08/2006

%%%%%%%%%%%%%%%%
% INITIALIZATION
%
intercept=ones(length(x),1);
geg=[x,y];
[n,m]=size(geg);
alfa=0.75;
hdefault=min(floor(2*floor((n+m+1)/2)-n+2*(n-floor((n+m+1)/2))*alfa),n);

if nargin < 3
    options.alpha=alfa;
    options.h=hdefault;
    options.ntrial=500;
    options.plots=1;
    options.classic=0;
    options.Hsets=[];
else
    default=struct('alpha',alfa,'h',hdefault,'ntrial',500,'plots',1,'classic',0,'Hsets',[]);
    list=fieldnames(default);
    options=default;
    IN=length(list);
    i=1;   
    counter=1;
    %   
    if nargin>2
        %
        %placing inputfields in array of strings
        %
        for j=1:nargin-2
            if rem(j,2)~=0
                chklist{i}=varargin{j};
                i=i+1;
            end
        end 
        dummy=sum(strcmp(chklist,'h')+2*strcmp(chklist,'alpha'));
        switch dummy
            case 0 %no input for alpha or h so take on the default values 
                options.alpha=0.75;
                options.h=floor(options.alpha*n);
            case 3
                error('Both input arguments alpha and h are provided. Only one is required.')
        end
        %
        %Checking which default parameters have to be changed
        % and keep them in the structure 'options'.
        %
        while counter<=IN 
            index=strmatch(list(counter,:),chklist,'exact');
            if ~isempty(index) %in case of similarity
                for j=1:nargin-2 %searching the index of the accompanying field
                    if rem(j,2)~=0 %fieldnames are placed on odd index
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
        if dummy==1
            options.alpha=options.h/n;
        elseif dummy==2
            options.h=floor(options.alpha*n);
        end
        Hsets = options.Hsets;
    end
end
%%%%%%%%
% MAIN %
[mcd_res, mcd_raw]=mcdcov(geg,'h',options.h,'plots',0,'Hsets',options.Hsets);
options.h = mcd_res.h;

%in case of an exact fit, the calculations stop at this point.
if ~isempty(mcd_res.plane)
    disp('Warning (mcdregres): The MCD covariance matrix is singular. See also mcdcov.m.')
    result=mcd_res;
    return
end

mcdreg.Hsubsets.Hopt = mcd_res.Hsubsets.Hopt;
mcdreg.Hsubsets.Hfreq = mcd_res.Hsubsets.Hfreq;

rewcovmcd=mcd_res.cov; %one-step reweighted covariance matrix
rewcenmcd=mcd_res.center; %one-step reweighted location
q=size(y,2);
p=size(x,2);

%initializing reweighted location and reweighted scatter (paragraph 5.1 of paper)
rewtmcdx=rewcenmcd(1:p)'; %column vectors!!!
rewtmcdy=rewcenmcd((p+1):m)';

rewsmcdx=rewcovmcd(1:p,1:p);
rewsmcdxy=rewcovmcd(1:p,(p+1):m);
rewsmcdyx=rewsmcdxy';
rewsmcdy=rewcovmcd((p+1):m,(p+1):m);

covyy=rewsmcdy;
covxx = rewsmcdx;
covxy = rewsmcdxy;
covyx = rewsmcdyx;
mcdreg.Sigma = [covxx, covxy ; covyx, covyy];
mcdreg.Mu = [rewtmcdx;rewtmcdy];

%reweighted beta and alpha (the columnvector [\beta^L ; \alpha^L])
rewbetamcd=[inv(rewsmcdx)*rewsmcdxy; (rewtmcdy-(rewsmcdyx*inv(rewsmcdx)*rewtmcdx))'];

%calculation of the reweighted weights based on 
%the residuals calculated with the reweighted beta-coefficients 
rewweights=zeros(n,1);
weights=zeros(n,1);
rewfitted=[x,intercept]*rewbetamcd(1:(p+1),:);
rewresid=y-rewfitted; %(r_i^L)
rewE=rewsmcdy-rewbetamcd(1:p,1:q)'*rewsmcdx*rewbetamcd(1:p,1:q); %reweighted scatter (\Sigma_eps^L)
for j=1:n
    if (sqrt(rewresid(j,1:q)*inv(rewE)*rewresid(j,1:q)')) <= sqrt(chi2inv(0.99,q))
        rewweights(j)=1;
    end
end

%regression reweighting part based on the reweighted observations. (paragraph 5.3 of paper)
rewclasscov=cov(geg(rewweights==1,:));
rewclasscenter=mean(geg(rewweights==1,:));

rewtmcdx=rewclasscenter(1:p)';
rewtmcdy=rewclasscenter((p+1):m)';

rewsmcdx=rewclasscov(1:p,1:p);
rewsmcdxy=rewclasscov(1:p,(p+1):m);
rewsmcdyx=rewsmcdxy';
rewsmcdy=rewclasscov((p+1):m,(p+1):m);

%regression reweighted coefficients beta and alpha ([\beta^{RL} \alpha^{RL}])
rewbetamcdrew=[inv(rewsmcdx)*rewsmcdxy; (rewtmcdy-(rewsmcdyx*inv(rewsmcdx)*rewtmcdx))'];
%regression reweighted scatter (\Sigma_eps^{RL}])
rewE2=rewsmcdy-rewbetamcdrew(1:p,1:q)'*rewsmcdx*rewbetamcdrew(1:p,1:q);
rewfittedrew=[x,intercept]*rewbetamcdrew(1:(p+1),:);
rewresidrew=y-rewfittedrew;

%%%%%%%%%%%%%%%%%
%OUTPUT STRUCTURE
%
mcdreg.covyy=rewsmcdy; %regression and location reweighted covariance matrix of responses
mcdreg.covxx = rewsmcdx;
mcdreg.covxy = rewsmcdxy;
mcdreg.covyx = rewsmcdyx;
mcdreg.Sigmarew = [mcdreg.covxx, mcdreg.covxy ; mcdreg.covyx, mcdreg.covyy];
mcdreg.Murew = [rewtmcdx;rewtmcdy];

mcdreg.x=x;
mcdreg.y=y;
mcdreg.coeffs=rewbetamcdrew; %regression and location reweighted coefficients
mcdreg.cov=rewE2; %scatter matrix based on location and regression reweighting
mcdreg.fitted=rewfittedrew; %estimated respons(es);
mcdreg.res=rewresidrew; %regression and location reweighted residuals

if(intercept)
   mcdreg.interc=1;
else
   mcdreg.interc=0;
end

% Robust distances in x-space = x-distances (rd) needed in diagnostic regression plot
if (-log(det(mcd_res.cov))/m) > 50
   mcdreg.rd='singularity';
else
   mcdreg.rd=sqrt(mahalanobis(mcdreg.x,rewcenmcd(1:p),'cov',rewcovmcd(1:p,1:p)))';
end

% Robust residual distances (resd) needed in diagnostic regression plot 
if q>1
   if (-log(det(mcd_res.cov))/m)>50
       mcdreg.resd='singularity';
       disp('Warning (mcdregres): A singularity ')
   else
      cen=zeros(q,1)';
      [nn,pp]=size(rewresidrew);
      mcdreg.resd=sqrt(mahalanobis(rewresidrew,cen,'cov',mcdreg.cov))'; %robust distances of residuals
  end
else
    mcdreg.covarRes=sqrt(mcdreg.cov);
    mcdreg.resd=rewresidrew(:,1)/mcdreg.covarRes; %standardized residuals 
end

% cutoff values 
mcdreg.cutoff.rd=sqrt(chi2inv(0.975,p)); 
mcdreg.cutoff.resd=sqrt(chi2inv(0.975,q));
mcdreg.flag=(abs(mcdreg.resd)<=mcdreg.cutoff.resd);

% robust multivariate Rsquared
mcdreg.weights=rewweights;
Yw=y(mcdreg.weights==1,:);
cYw=mcenter(Yw);
res=rewresidrew(mcdreg.weights==1,:);
mcdreg.rsquared=1-(det(res'*res)/det(cYw'*cYw));
mcdreg.class='MCDREG';


if options.classic
    mcdreg.classic=mlr(x,y,'plots',0);
else
    mcdreg.classic=0;
end

result=struct('slope',{mcdreg.coeffs(1:p,:)}, 'int',{mcdreg.coeffs(p+1,:)}, ...
    'fitted',{mcdreg.fitted},'res',{mcdreg.res},'cov',{mcdreg.cov},'rsquared',{mcdreg.rsquared},...
    'h',{options.h},'Hsubsets',{mcdreg.Hsubsets},'rd', {mcdreg.rd},'resd',{mcdreg.resd},'cutoff',{mcdreg.cutoff},...
    'weights',{mcdreg.weights},'flag',{mcdreg.flag'},'class',{mcdreg.class},'classic',{mcdreg.classic},...
    'Mu',{mcdreg.Mu},'Sigma',{mcdreg.Sigma},'Murew',{mcdreg.Murew},'Sigmarew',{mcdreg.Sigmarew});

try
    if options.plots & options.classic
        makeplot(result,'classic',1)
    elseif options.plots
        makeplot(result)    
    end
catch %output must be given even if plots are interrupted 
    %> delete(gcf) to get rid of the menu 
end