function result=cpcr(x,y,varargin)

%CPCR performs a classical principal components regression.
% First, classical PCA is applied to the predictor variables x (see cpca.m) and 
% k components are retained. Then a multiple linear regression method (see mlr.m)
% is performed of the response variable y on the k principal components. 
%
% I/O: result=cpcr(x,y,'k',2);
%
% Required input arguments:
%      x : Data matrix of the explanatory variables
%          (n observations in rows, p variables in columns)
%      y : Data matrix of the response variables
%          (n observations in rows, q variables in columns)
%
% Optional input argument: 
%      k : Number of principal components to compute. If k is missing, 
%          a scree plot is drawn which allows to select
%          the number of principal components.
%  plots : If equal to one (default), a menu is shown which allows to draw several plots,
%          such as a score outlier map and a regression outlier map. 
%          If 'plots' is equal to zero, all plots are suppressed.
%          See also makeplot.m
%     
% The output is a structure containing 
%
%   result.slope      : Classical slope
%   result.int        : Classical intercept
%   result.fitted     : Classcial prediction vector
%   result.res        : Classical residuals
%   result.sigma      : Estimated variance-covariance matrix of the residuals 
%   result.rsquared   : R-squared value
%   result.k          : Number of components used in the regression
%   result.sd         : Classical score distances
%   result.od         : Classical orthogonal distances
%   result.resd       : Residual distances (when there are several response variables).
%                       If univariate regression is performed, it contains the standardized residuals.
%   result.cutoff     : Cutoff values for the score, orthogonal and residual distances.
%   result.flag      : The observations whose orthogonal distance is larger than result.cutoff.od
%                      (orthogonal outliers => result.flag.od) and/or whose residual distance is 
%                      larger than result.cutoff.resd (bad leverage points/vertical outliers => result.flag.resd)
%                      can be considered as outliers and receive a flag equal to zero (=> result.flag.all).
%                      The regular observations, including the good leverage points, receive a flag 1.
%   result.class      : 'CPCR'
%   result.cpca       : Full output of the classical PCA part (see cpca.m)
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by S. Verboven
% Last Update: 05/04/2003 

default=struct('plots',1,'k',0);
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


%classical PCA
if options.k==0
    out.cpca=cpca(x,'plots',0);
else
    out.cpca=cpca(x,'plots',0,'k',options.k);
end
%model with intercept 
T=[out.cpca.T(:,1:out.cpca.k) ones(size(out.cpca.T,1),1)];
%Multivariate linear regression
out.a=inv(T'*T)*T'*y;
out.fitted=T*out.a;
len=size(out.a,1);
p=size(T,2);
q=size(y,2);
geg=[T,y];
[n,m]=size(geg);
S=cov(geg);
Sx=S(1:p-1,1:p-1);
Sxy=S(1:p-1,p+1:m);
Syx=Sxy';
Sy=S(p+1:m,p+1:m);
Se=Sy-out.a(1:p-1,1:q)'*Sx*out.a(1:p-1,1:q); %variance of errors
%regression coefficients slope and intercept ([\beta \alpha])
out.coeffs=[out.cpca.P(:,1:out.cpca.k)*out.a(1:len-1,:); out.a(len,:)-out.cpca.M*out.cpca.P(:,1:out.cpca.k)*out.a(1:len-1,:)];    %%coefficients in the original space;
out.slope=out.cpca.P(:,1:out.cpca.k)*out.a(1:len-1,:);
out.int=out.a(len,:)-out.cpca.M*out.cpca.P(:,1:out.cpca.k)*out.a(1:len-1,:);
out.res=y-out.fitted;
out.sigma=Se;
if q==1
    out.stdres=out.res./sqrt(diag(out.sigma));
end
out.k=out.cpca.k;
STTm=sum((y-repmat(mean(y),length(y),1)).^2);
SSE=sum(out.res.^2);
out.rsquared=1-SSE/STTm;
out.class='CPCR';

%calculating residual distances
if q>1
    if (-log(det(Se)/(m-1)))>50
        out.resd='singularity';
    else
       cen=zeros(q,1);
       out.resd=sqrt(mahalanobis(out.res,cen','cov',Se))';
   end
else % q==1
    out.resd=out.stdres; %standardized residuals 
end
out.cutoff=out.cpca.cutoff;
out.cutoff.resd=sqrt(chi2inv(0.975,size(y,2)));
%computing flags
out.flag.od=out.cpca.flag.od;
out.flag.sd=out.cpca.flag.sd;
out.flag.resd=abs(out.resd)<=out.cutoff.resd;
out.flag.all=(out.flag.od & out.flag.resd);

result=struct( 'slope',{out.slope}, 'int',{out.int},'fitted',{out.fitted},'res',{out.res},...
    'cov',{out.sigma},'rsquared',{out.rsquared},...
    'k',{out.k},'sd', {out.cpca.sd},'od',{out.cpca.od},'resd',{out.resd},...
    'cutoff',{out.cutoff},'flag',{out.flag},'class',{out.class},...
    'cpca',{out.cpca});

try
    if options.plots
        makeplot(result)
    end
catch %output must be given even if plots are interrupted 
    %> delete(gcf) to get rid of the menu 
    end
    