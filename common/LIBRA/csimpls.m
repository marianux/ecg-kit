function result=csimpls(x,y,varargin)

%CSIMPLS performs Partial Least Squares regression using the SIMPLS 
% algorithm of de Jong (1993).
%
% Reference: 
%    de Jong, S. (1993),
%    "SIMPLS: an alternative approach to Partial Least Squares Regression", 
%    Chemometrics and Intelligent Laboratory Systems, 18, 251-263.
%
% Required input arguments:
%            x : Data matrix of the explanatory variables
%                (n observations in rows, p variables in columns)
%            y : Data matrix of the response variables
%                (n observations in rows, q variables in columns)
%     
% Optional input arguments:
%            k : Number of components to be used.  
%                (default = min([9,rank([x,y]),floor(n/2),p])). 
%        plots : If equal to one, a menu is shown which allows to draw several plots,
%                such as a score outlier map and a regression outlier map. (default)
%                If 'plots' is equal to zero, all plots are suppressed.
%                See also makeplot.m
%
% I/O: result=csimpls(x,y,'k',10);
%
% The output of CSIMPLS is a structure containing:
%
%   result.slope     : Classical slope estimate
%   result.int       : Classical intercept estimate
%   result.fitted    : Classical prediction vector
%   result.res       : Classical residuals
%   result.cov       : Estimated variance-covariance matrix of the residuals 
%   result.M         : Classical center of the matrix [X;Y]
%   result.T         : Classical scores
%   result.weights.r : Classical simpls weights
%   result.weights.p : Classical simpls weights
%   result.Tcov      : Classical covariance matrix of the scores
%   result.k         : Number of components used in the regression
%   result.sd        : Classical score distances
%   result.od        : Classical orthogonal distances
%   result.resd        : Residual distances (when there are several response variables).
%                      If univariate regression is performed, it contains the standardized residuals.
%   result.cutoff    : Cutoff values for the score, orthogonal and residual distances.
%   result.flag      : The observations whose orthogonal distance is larger than result.cutoff.od
%                      (orthogonal outliers => result.flag.od) and/or whose residual distance is 
%                      larger than result.cutoff.resd (bad leverage points/vertical outliers => result.flag.resd)
%                      can be considered as outliers and receive a flag equal to zero (=> result.flag.all).
%                      The regular observations, including the good leverage points, receive a flag 1.
%   result.class     : 'CSIMPLS'
%
% This function is part of LIBRA: the Matlab library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
%Written by Karlien Vanden Branden
%Last Update: 05/04/2003 
%Last revision: 20/05/2008

[n,p1]=size(x);
[n2,q1]=size(y);
if n~=n2
    error('The response variables and the predictor variables have a different number of observations.')
end
kmin=min([9,rank(x),floor(n/2),p1]);
default=struct('plots',1,'k',kmin);
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

[xcentr,cx]=mcenter(x);
[ycentr,cy]=mcenter(y);
out.M=[cx cy];%[xcentr ycentr];
sigmaxy=xcentr'*ycentr;
for i=1:options.k
    sigmayx=sigmaxy';
    if q1>p1        
        [RR,LL]=eig(sigmaxy*sigmayx);  
        [LL,I]=greatsort(diag(LL));
        rr=RR(:,I(1));
        qq=sigmayx*rr; 
        qq=qq/norm(qq);
    else
        if q1==1
            qq = 1;
            rr = sigmaxy;
        else
            [QQ,LL]=eig(sigmayx*sigmaxy);
            [LL,I]=greatsort(diag(LL));
            qq=QQ(:,I(1));
            rr=sigmaxy*qq;
        end
    end
    tt=xcentr*rr;
    nttc=norm(tt);
    rr=rr/nttc;
    tt=tt/nttc;
    qq=ycentr'*tt;
    uu=ycentr*qq;
    pp=xcentr'*tt;
    vv=pp;
    if i>1 												
        vv = vv -v*(v'*pp);
    end
    vv = vv/norm(vv);
    sigmaxy = sigmaxy - vv*(vv'*sigmaxy);
    v(:,i)=vv;
    q(:,i)=qq;
    t(:,i)=tt;
    u(:,i)=uu;
    p(:,i)=pp;
    r(:,i)=rr;
end

%Second Stage : Classical Regression and transformation

b=r*q';
int=cy-cx*b;

%classical output
out.T=t;
out.weights.p=p;
out.weights.r=r;
out.coef=[b; int];
out.b=b;
out.int=int;
out.yhat=x*b+repmat(int,n,1);
out.res=y-out.yhat;
out.covar=cov(out.res);
if q1==1
    out.stdres=out.res./sqrt(out.covar);
end
out.k=options.k;
STTm=sum((y-repmat(mean(y),length(y),1)).^2);
SSE=sum(out.res.^2);
out.rsquare=1-SSE/STTm;
out.class='CSIMPLS';

%calculating classical distances in x-space
out.Tcov=cov(t);
out.sd=sqrt(mahalanobis(t,zeros(1,out.k),'cov',out.Tcov))';
out.cutoff.sd=sqrt(chi2inv(0.975,out.k));

%calculating classical orthogonal distances
xt=t*p';
tempo=xcentr-xt;
for i=1:n
    out.od(i,1)=norm(tempo(i,:));
end
r=rank(x);
if out.k~=r
    m=mean(out.od.^(2/3));
    s=sqrt(var(out.od.^(2/3)));
    out.cutoff.od = sqrt(norminv(0.975,m,s)^3); 
else
    out.cutoff.od=0;
end

%calculating residual distances
if q1>1
    if (-log(det(out.covar)/(p1+q1-1)))>50
        out.resd='singularity';
    else
       cen=zeros(q1,1);
       out.resd=sqrt(mahalanobis(out.res,cen','cov',out.covar))';
   end
else    %here q==1
    out.resd=out.stdres; %standardized residuals
end
out.cutoff.resd=sqrt(chi2inv(0.975,q1));

%Computing flags
out.flag.od=out.od<=out.cutoff.od;
out.flag.resd=abs(out.resd)<=out.cutoff.resd;
out.flag.all=(out.flag.od & out.flag.resd);


result=struct('slope',{out.b}, 'int',{out.int},'fitted',{out.yhat},'res',{out.res}, 'cov',{out.covar},...
    'M', {out.M},'T',{out.T} ,'weights', {out.weights},'Tcov',{out.Tcov},'k',{out.k},'sd', {out.sd},'od',{out.od},'resd',{out.resd},...
    'cutoff',{out.cutoff},'flag',{out.flag},'class',{out.class});

try
    if options.plots
        makeplot(result)
    end
catch %output must be given even if plots are interrupted 
end






