function result = mlr(x,y,varargin)

%MLR is the classical least squares estimator for multivariate multiple
% linear regression. It can handle both one or several ('multiple')
% predictor variables and one or several ('multivariate') response variables. 
% If the regression model contains an intercept, the x-matrix may not contain
% a column with ones. 
% If there is only one response variable, the function ols.m is called. 
%
% See for example: 
%   R.A. Johnson and D.W. Wichern, 
%   "Applied multivariate statistical analysis (Fifth Edition)",
%   Prentice Hall, chapter 7.
%
% Required input arguments:
%         x : Data matrix of the explanatory variables
%             (n observations in rows, p variables in columns)
%         y : Data matrix of the response variables
%             (n observations in rows, q variables in columns)
%
% Optional input arguments: 
% intercept : logical flag: if 1, a model with constant term will be 
%             fitted; if 0, no constant term will be included. (default: 1)
%     plots : If equal to one, a menu is shown which allows to draw several plots,
%             such as residual plots and a regression outlier map. (default)
%             If 'plots' is equal to zero, all plots are suppressed.
%             See also makeplot.m
%
% I/O: result=mlr(x,y,'plots',0);
%   The user should only give the input arguments that have to change their default value.
% 
% The output is a structure containing:
%
%   result.slope     : Slope estimate
%   result.int       : Intercept estimate
%   result.fitted    : Fitted values
%   result.res       : Residuals
%   result.stdres    : Standardized residuals 
%   result.cov       : Estimated variance-covariance matrix of the residuals
%   result.rsquared  : R-squared value
%   result.md        : Score distances (Mahalanobis distances in x-space)
%   result.resd      : Residual distances (when there are several response variables).
%                      If univariate regression is performed, it contains the standardized residuals.
%   result.cutoff    : Cutoff values for the score and residual distances
%   result.flag      : The observations whose residual distance is larger than result.cutoff.resd
%                      receive a flag equal to zero. The other observations receive a flag 1.
%   result.class     : 'MLR' (when q > 1) or 'LS' (when q = 1)
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Nele Smets on 13/12/2003
% Last update: 05/04/2004
%

q=size(y,2);
p=size(x,2);
geg=[x,y];
[n,m]=size(geg);
intercept=ones(n,1);
default=struct('plots',1,'intercept',1);
list=fieldnames(default);
options=default;
IN=length(list);
i=1;   
counter=1;
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
    % Checking which default parameters have to be changed
    % and keep them in the structure 'options'.
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
end
%%%%%%%%%Main part%%%%%%%
if q==1 
    result=ols(x,y,'intercept', options.intercept,'plots',options.plots);
else
    cmean=mean(geg); 
    ccovar=cov(geg);
    meanx=cmean(1:p)';
    meany=cmean((p+1):m)';
    cy=mcenter(y);   
    covarx=ccovar(1:p,1:p);
    covary=ccovar((p+1):m,(p+1):m);
    covarxy=ccovar(1:p,(p+1):m);
    covaryx=covarxy';
    res.betas=[inv(covarx)*covarxy; (meany-(covaryx*inv(covarx)*meanx))'];
    res.fitted=[x intercept]*res.betas;
    res.residuals=y-res.fitted;
    res.cov=covary-res.betas(1:p,1:q)'*covarx*res.betas(1:p,1:q);
    res.stdresid= res.residuals./repmat(diag(res.cov)',n,1);
    res.rsquared= 1-(det(res.residuals'*res.residuals)/det(cy'*cy));
    res.class='MLR';
    
    % x-distances (md) needed in diagnostic regression plot
    if (-log(det(ccovar))/m) > 50
        res.md='singularity';
        res.resd='singularity';
    else
        res.md=sqrt(mahalanobis(x,cmean(1:p),'cov',covarx))';
        cen=zeros(q,1)';
        res.resd=sqrt(mahalanobis(res.residuals,cen,'cov',res.cov))'; 
    end
    
    %cutoff values 
    res.cutoff.md=sqrt(chi2inv(0.975,p)); 
    res.cutoff.resd=sqrt(chi2inv(0.975,q));
    res.flag=(abs(res.resd)<=res.cutoff.resd);
    result=struct('slope',{res.betas(1:p,:)},'int',{res.betas(p+1,:)}, ...
        'fitted',{res.fitted},'res',{res.residuals},'stdres',{res.stdresid},...
        'cov',{res.cov},'rsquared',{res.rsquared},'md',{res.md},'resd',{res.resd},...
        'cutoff',{res.cutoff},'flag',{res.flag},'class',{res.class});
    if ~intercept
        result=setfield(result, 'int', 0)
    end
    try
        if options.plots==1
            makeplot(result)
        end
    catch %output must be given even if plots are interrupted 
        %> delete(gcf) to get rid of the menu 
    end
end




