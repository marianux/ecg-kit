function result = ols(x,y,varargin)

%OLS is the classical least squares estimator for multiple
% linear regression. It can handle both one or several predictor variables,
% and one response variable. 
% If there are several response variables, the function mlr.m should be used.
%
% Required input arguments:
%           x : Data matrix of the explanatory variables
%               (n observations in rows, p variables in columns)
%               Missing values (NaN's) and infinite values (Inf's) are allowed, since observations (rows) 
%               with missing or infinite values will automatically be excluded from the computations.
%           y : Data vector with the response variable
%               Missing values (NaN's) and infinite values (Inf's) are allowed, since observations (rows) 
%               with missing or infinite values will automatically be excluded from the computations.
%
% Optional input arguments:
%    intercept : logical flag: if 1, a model with constant term will be 
%                fitted; if 0, no constant term will be included. (default: 1)
%        plots : if 0, the plots are supressed (default:0)
%
% I/O:  result=ols(x,y,'plots',0,'intercept',0)
%  The user should only give the input arguments that have to change their default value.
%  The name of the input arguments needs to be followed by their value.
%  The order of the input arguments is of no importance.
%
% The output is a structure containing:
%
%   result.slope     : Slope estimate
%   result.int       : Intercept estimate (if no intercept is included, it equals zero)
%   result.fitted    : Fitted values
%   result.res       : Residuals
%   result.scale     : Scale estimate of the residuals
%   result.rsquared  : R-squared value
%   result.md        : Mahalanobis distances in x-space
%   result.resd      : Residual distances (which are equal to the standardized residuals)
%   result.cutoff    : Cutoff values for the score distances, and for the standardized residuals 
%   result.flag      : The observations whose absolute standardized residual is larger than result.cutoff.resd
%                      receive a flag equal to zero. The other observations receive a flag 1.
%   result.X         : If x is univariate, data matrix without missing or infinite values.
%   result.y         : If x is univariate, response vector without missing or infinite values.
%   result.class     : 'LS'
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Nele Smets on 06/12/2003
% Last update on 05/04/2004

if nargin<2
    error('there is a missing input argument')
end
default=struct('intercept',1,'plots',1);
list=fieldnames(default);
options=default;
IN=length(list);
i=1;
counter=1;
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
end
intercept=options.intercept;
plots=options.plots;


[n,p]=size(x);
na.x=~isfinite(x*ones(p,1));
na.y=~isfinite(y);
if size(na.x,1)~=size(na.y,1)
    error('Number of observations in x and y are not equal.');
end
ok=~(na.x|na.y);
x=x(ok,:);
y=y(ok,:);
n=length(y);
X=x;
if intercept
    x=cat(2,x,ones(n,1));
    p=p+1;
end

[coeff,bint,res] = regress(y,x);
fitted=x*coeff;
scale=sqrt(1/(n-p)*sum(res.^2));
stdres=res/scale;
md=sqrt(mahalanobis(X,mean(X),'cov',cov(X)))';
SSE=sum((y-fitted).^2);
if intercept
    SST=sum((y-mean(y)).^2);
    cutoff.md=sqrt(chi2inv(0.975,p-1));
else 
    SST=sum(y.^2);
    cutoff.md=sqrt(chi2inv(0.975,p));
end
cutoff.resd=sqrt(chi2inv(0.975,1));
rsquared=1-SSE/SST;
flags=(abs(stdres)<=cutoff.resd);

result=struct('slope',{coeff(1:p)},'int',{0},'fitted',{fitted},'res',{res},'scale',{scale},'rsquared',{rsquared},...
    'md',md,'resd', {stdres},'cutoff',cutoff,'flag',{flags},'class',{'LS'},'X',{X},'y',{y});

if intercept
    result=setfield(result,'slope',coeff(1:p-1));
    result=setfield(result,'int', coeff(p));
end
if size(X,2)~=1
    result=rmfield(result,{'X','y'});
end
try
    if plots
        makeplot(result)
    end
catch %output must be given even if plots are interrupted 
    %> delete(gcf) to get rid of the menu 
end
