function result=rsimpls(x,y,varargin)

%RSIMPLS is a 'Robust method for Partial Least Squares Regression based on the
% SIMPLS algorithm'. It can be applied to both low and high-dimensional predictor variables x
% and to one or multiple response variables y. It is resistant to outliers in the data.
% The RSIMPLS algorithm is built on two main stages. First, a matrix of scores is derived 
% based on a robust covariance criterion (see robpca.m),
% and secondly a robust regression is performed based on the results from ROBPCA. 
% 
% The RSIMPLS method is described in: 
%    Hubert, M., and Vanden Branden, K. (2003),
%    "Robust Methods for Partial Least Squares Regression",
%    Journal of Chemometrics, 17, 537-549.
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
%            k : Number of components to be used. 
%                (default = min(rank([x,y]),kmax)). If k is not specified,
%                it can be selected using the option 'rmsecv'. 
%         kmax : Maximal number of components to be used (default = 9).
%                If k is provided, kmax does not need to be specified, unless k is larger
%                than 9.                 
%        alpha : (1-alpha) measures the fraction of outliers the algorithm should 
%                resist. Any value between 0.5 and 1 may be specified. (default = 0.75)
%            h : (n-h+1) measures the number of observations the algorithm should 
%                resist. Any value between n/2 and n may be specified. (default = 0.75*n)
%                Alpha and h may not both be specified.
%       rmsecv : If equal to zero and k is not specified, a robust R-squared curve
%                is plotted and an optimal k value can be chosen. (default)
%                If equal to one and k is not specified, a robust component selection-curve is plotted and
%                an optimal k value can be chosen. This curve computes a combination of
%                the cross-validated root mean squared error and the robust residual 
%                sum of squares for k=1 to kmax. Rmsecv and k may not both
%                be specified.
%        rmsep : If equal to one, the robust RMSEP-value (root mean squared error of
%                prediction) for the model with k components (default = 0). This value
%                is automatically given if rmsecv = 1.
%        plots : If equal to one, a menu is shown which allows to draw several plots,
%                such as a robust score outlier map and a regression
%                outlier map are drawn. (default)
%                If the input argument 'classic' is equal to one, the classical
%                diagnostic plots are drawn as well.
%                If 'plots' is equal to zero, all plots are suppressed.
%                See also makeplot.m
%        labsd : The 'labsd' observations with largest score distance are
%                labeled on the outlier map (default = 3)
%        labod : The 'labod' observations with largest orthogonal distance are
%                labeled on the outlier map (default = 3)    
%      labresd : The 'labresd' observations with largest residual distance are
%                labeled on the outlier map (default = 3)   
%      classic : If equal to one, the classical SIMPLS analysis will be performed as well
%                (see also csimpls.m). (default = 0)
%
% Options for advanced users:
%           kr : Total number of components used by the ROBPCA method.
%                We advise to use kr=k+q. (default) 
%        kmaxr : Maximal number of components used by the ROBPCA method. 
%                default = min(kmax+q,rank([x,y]))
%  plotsrobpca : If equal to one, a robust score outlier map from ROBPCA is drawn. 
%                If the input argument 'classic' is equal to one, the classical
%                outlier map is drawn as well.
%                If 'plotsrobpca' is equal to zero, all plots are suppressed. (default)
%           st : Indicates the current stage of the algorithm for cross-validation (RMSECV)  
%                default = 0 -> performs all the stages of the algorithm and computes all the 
%                parameters for the full model. 
%                st=1, robpca is still performed, but when 
%                st=2, the algorithm proceeds based on the previous knownledge of the output from robpca. 
%          out : Is an empty structure, but is constructed while performing the cross-validation.
%                (see rrmse.m)
%
% I/O: result=rsimpls(x,y,'k',k,'kmax',10,'alpha',0.75,'h',h,'rmsecv',0,'rmsep',0,...
%       'plots',1,'labsd',3,'labod',3,'labresd',3,'classic',1,'kr',kr,...
%       'kmaxr',kmaxr,'plotsrobpca',0,'st',0,'out',[]);
%  The user should only give the input arguments that have to change their default value.
%  The name of the input arguments needs to be followed by their value.
%  The order of the input arguments is of no importance.
%
% Examples:
%   rsimpls(x,y,'k',5,'plots',1);
%   rsimpls(x,y,'classic',1,'rmsecv',1);
%
% The output of RSIMPLS is a structure containing:
%
%   result.slope     : Robust slope estimate
%   result.int       : Robust intercept estimate
%   result.fitted    : Robust fitted values
%   result.res       : Robust residuals
%   result.cov       : Estimated variance-covariance matrix of the residuals 
%   result.T         : Robust scores
%   result.weights.r : Robust simpls weights
%   result.weights.p : Robust simpls weights
%   result.Tcenter   : Robust center of the scores
%   result.Tcov      : Robust covariance matrix of the scores 
%   result.rsquared  : Robust R-squared value for the optimal k
%   result.rcs       : Robust Component Selection Criterion:
%                      This is a matrix with kmax columns. The first row contains the approximate
%                      R-squared values (for k=1,...,kmax) and the second row the square root of the weighted
%                      residuals sum of squares. This is equal to the RCS-value with 
%                      gamma = 0 (see Engelen and Hubert (2005) for the definition). 
%                      If the input argument rmsecv = 1, the third row contains the RCS-value for
%                      gamma = 0.5 and the fourth for gamma = 1. The last one is equal to the robust 
%                      cross-validated RMSE values. 
%                      Note that all the entries in this matrix depend on the choice of kmax. 
%   result.rmsep     : Robust RMSEP value 
%   result.k         : Number of components used in the regression
%   result.h         : The quantile h used throughout the algorithm
%   result.sd        : Robust score distances
%   result.od        : Robust orthogonal distances
%   result.resd      : Residual distances (when there are several response variables).
%                      If univariate regression is performed, it contains the standardized residuals.
%   result.cutoff    : Cutoff values for the score (result.cutoff.sd), orthogonal 
%                      (result.cutoff.od) and residual distances (result.cutoff.resd).
%                      We use 0.975 quantiles of the chi-squared distribution.
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
%   result.class     : 'RSIMPLS'
%   result.classic   : If the inputargument 'classic' is equal to 1, this structure
%                      contains results of the classical SIMPLS analysis. (see also csimpls.m)
%   results.robpca   : The results of robpca on [X,Y]. 
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Karlien Vanden Branden
% Version date: 07/04/2004 
% Last update: 04/08/2006

%
%initialization with defaults
%
if rem(nargin,2)~=0
    error('Number of input arguments must be even!');
end
[n,p1]=size(x);
[n2,q1]=size(y);
z=[x,y];
rx=rank(x);
rz=rank(z);
if n~=n2
    error('The response variables and the predictor variables have a different number of observations.')
end
niter=100;counter=1;mcd=0;alfa=0.75;
kmax=min([9,rx,floor(n/2),p1]);
kmaxr=min([kmax+q1,rz]);
h=floor(2*floor((n+kmaxr+1)/2)-n+2*(n-floor((n+kmaxr+1)/2))*alfa);
labsd=3;labod=3;labresd=3;
plotsrobpca=0;plots=1;k=0;
kr=k+q1;
st=0;rmsecv=0;
out=[];classic=0;rmsep=0;rmsep_value=nan;rmsecv_value = nan;rsquared_value = nan;rss_value = nan;
default=struct('alpha',alfa,'h',h,'labsd',labsd,'labod',labod,'labresd',labresd,'k',k,'kr',kr,...
    'plotsrobpca',plotsrobpca,'plots',plots,'kmax',kmax,'st',st,...
    'out',out,'rmsecv',rmsecv,'classic',classic,'rmsep',rmsep,...
    'kmaxr',kmaxr,'rmsep_value',rmsep_value,'rmsecv_value',rmsecv_value,'rsquared_value',...
    rsquared_value,'rss_value',rss_value);
list=fieldnames(default);
options=default;
IN=length(list);
i=1;
%
if nargin>2
    %
    %placing inputfields in array of strings
    %
    for j=1:nargin-3
        if rem(j,2)~=0
            chklist{i}=varargin{j};
            i=i+1;
        end
    end 
    dummy=sum(strcmp(chklist,'h')+2*strcmp(chklist,'alpha'));
    switch dummy
        case 0 %Take on default values  
            options.alpha=alfa;%0.75;
            options.h=h;
        case 3
            error('Both input arguments alpha and h are provided. Only one is required.')
    end
    %
    %Checking which default parameters have to be changed
    % and keep them in the structure 'options'.
    %
    while counter<=IN 
        index=strmatch(list(counter,:),chklist,'exact');%contains the users input one by one
        if ~isempty(index) %in case of similarity
            for j=1:nargin-3 %searching the index of the accompanying field
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
    options.h=floor(options.h);
    options.kmax=floor(options.kmax);
    options.k=floor(options.k);
    options.kmaxr=floor(options.kmaxr);
    options.kr=floor(options.kr);
    kmax=min([options.kmax,floor(n/2),rz,p1]);
    kmaxr=max([min([options.kmaxr,rz]),kmax+q1]);
    k=min(options.k,kmax);
    while k<0
        k=input(['The number of components can not be negative.\n'...
        'How many principal components would you like to retain?\n']);
    end
    if any(strcmp(chklist,'kr'))
        kr=floor(max([options.kr,k+q1]));
    else
        kr=k+q1;
    end
    if dummy==1 %checking inputvariable h
        if options.h-floor(options.h)~=0
            mess=sprintf('Attention (rsimpls.m): h must be an integer. \n');
            disp(mess)
        end
        if kr==0
            if options.h<floor((n+kmaxr+1)/2) 
                options.h=floor((n+kmaxr+1)/2);
                mess=sprintf(['Attention (rsimpls.m): h should be larger than (n+kmaxr+1)/2.\n',...
                        'It is set to its minimum value',num2str(options.h)]);
                disp(mess)
            end
            options.alpha=options.h/n;
        else
            if options.h<floor((n+kr+1)/2) 
                options.h=floor((n+kr+1)/2);
                mess=sprintf(['Attention (rsimpls.m): h should be larger than (n+kr+1)/2.\n',...
                        'It is set to its minimum value',num2str(options.h)]);
                disp(mess)
            end
            options.alpha=options.h/n;
        end
        if options.h>n
            options.alpha=0.75;
            if kr==0
                options.h=floor(2*floor((n+kmaxr+1)/2)-n+2*(n-floor((n+kmaxr+1)/2))*options.alpha);
            else
                options.h=floor(2*floor((n+kr+1)/2)-n+2*(n-floor((n+kr+1)/2))*options.alpha);
            end    
            mess=sprintf(['Attention (rsimpls.m): h should be smaller than n. \n',...
                    'It is set to its default value ',num2str(options.h)]);
            disp(mess)
        end
    elseif dummy==2
        if options.alpha < 0.5
            options.alpha=0.5;
            mess=sprintf(['Attention (rsimpls.m): Alpha should be larger than 0.5. \n',...
                    'It is set to 0.5.']);
            disp(mess)
        end
        if options.alpha > 1
            options.alpha=0.75;
            mess=sprintf(['Attention (rsimpls.m): Alpha should be smaller than 1.\n',...
                    'It is set to 0.75.']);
            disp(mess)
        end
        if kr==0
            options.h=floor(2*floor((n+kmaxr+1)/2)-n+2*(n-floor((n+kmaxr+1)/2))*options.alpha);
        else
            options.h=floor(2*floor((n+kr+1)/2)-n+2*(n-floor((n+kr+1)/2))*options.alpha);
        end  
    end
    h=options.h;alfa=options.alpha;labsd=max(0,min(floor(options.labsd),n));
    dummyh = strcmp(chklist,'h');
    dummykmax = strcmp(chklist,'kmax');
    if all(dummyh == 0) && any(dummykmax)
        h = floor(2*floor((n+kmax+1)/2)-n+2*(n-floor((n+kmax+1)/2))*alfa);
    end
    labod=max(0,min(floor(options.labod),n));labresd=max(0,min(floor(options.labresd),n));
    plotsrobpca=options.plotsrobpca;plots=options.plots;
    st=options.st;
    out=options.out;
    rmsecv=options.rmsecv;
    classic=options.classic;
    rmsep=options.rmsep;
    rmsep_value = options.rmsep_value;
    rmsecv_value = options.rmsecv_value;
    rmsecv = options.rmsecv;
    rsquared_value = options.rsquared_value;
    rss_value = options.rss_value;
end

if q1==1 && k>=(h-2)
    mess=sprintf(['Attention (rsimpls.m): The number of components, k = ',num2str(k),...
            '\n is larger than our recommended maximum value of k = ',num2str(h-2)-1,'.']);
    disp(mess)
elseif q1>1 && k>=((h/q1)-(q1/2)-0.5)
    mess=sprintf(['Attention (rsimpls.m): The number of components, k = ',num2str(k),...
            '\n is larger than our recommended maximum value of k = ',num2str(floor((h/q1)-(q1/2)-1.5)),'.']);
    disp(mess)
end

%
%MAIN PART
%
% selection of number of components
if k == 0
    if rmsecv == 0
        [R2,final]=rsquared(x,y,kmax,'RSIMPLS',options.h);
        rss = final.rss;
        k = final.k;
        result=rsimpls(x,y,'k',k,'kr',k+q1,'h',h,'rmsecv',0,'rmsep',options.rmsep,'plots',0,'classic',classic,'rsquared_value',R2,'rss_value',rss); 
    else
        out=rrmse(x,y,h,kmax,'RSIMPLS',1);
        R2 = out.R2;
        rss = out.rss;
        k=out.k;
        rmsecv_value = out.rmsecv;
        pred = rrmse(x,y,h,kmax,'RSIMPLS',0,k,out.weight,out.res);
        result=rsimpls(x,y,'k',k,'kr',k+q1,'h',h,'rmsecv',0,'rmsep',0,'plots',0,'classic',classic,'rmsep_value',pred.rmsep,'rmsecv_value',rmsecv_value,'rsquared_value',R2,'rss_value',rss); 
    end
else
    if rmsecv
        error(['Both RMSECV and k were given.', ...
            'Please rerun your analysis with one of these inputs.  (see help file)'])
    end
%First stage: Obtain the scores T by first performing ROBPCA on z: 
    if st<=1
        out.robpca=robpca(z,'k',kr,'h',h,'plots',plotsrobpca,'kmax',kmaxr,'classic',classic,'mcd',mcd);
        out.h=h;
        out.centerz=out.robpca.M;
        out.sigmaxy=out.robpca.P(1:p1,:)*diag(out.robpca.L)*out.robpca.P(p1+1:p1+q1,:)';
        out.sigmax=out.robpca.P(1:p1,:)*diag(out.robpca.L)*out.robpca.P(1:p1,:)';
        out.xcentr=x-repmat(out.centerz(1:p1),n,1);
        out.ycentr=y-repmat(out.centerz(p1+1:p1+q1),n,1);
        out.weights2=out.robpca.flag.all; 
    end
    if st
        i=k;
    else 
        i=1;
    end
    while i<=k
        out.sigmayx=out.sigmaxy';
        if q1>p1        
            [RR,LL]=eig(out.sigmaxy*out.sigmayx); 
            [LL,I]=greatsort(diag(LL));
            rr=RR(:,I(1));
            qq=out.sigmayx*rr;  						                      
            qq=qq/norm(qq); 
        else
            [QQ,LL]=eig(out.sigmayx*out.sigmaxy);	
            [LL,I]=greatsort(diag(LL));
            qq=QQ(:,I(1));
            rr=out.sigmaxy*qq;
            rr=rr/norm(rr); 
        end
        tt=out.xcentr*rr;
        uu=out.ycentr*qq;	
        pp=out.sigmax*rr/(rr'*out.sigmax*rr);
        vv=pp;
        if i>1 												
            vv=vv-out.v*(out.v'*pp);
        end
        if vv'*vv==0
            error('The number of components is too large')
        end
        vv=vv./norm(vv);									
        out.sigmaxy=out.sigmaxy-vv*(vv'*out.sigmaxy); 
        out.v(:,i)=vv;
        out.q(:,i)=qq;
        out.t(:,i)=tt;
        out.u(:,i)=uu;
        out.p(:,i)=pp;
        out.r(:,i)=rr;
        i=i+1;
    end
    
    %Second Stage : Robust ROBPCA-regression
    robpcareg=robpcaregres(out.t,y,out.weights2);
    breg=robpcareg.coeffs(1:k,:);
    Yhat=out.t*breg+repmat(robpcareg.coeffs(k+1,:),n,1);
    b=out.r*breg;  
    int=robpcareg.coeffs(k+1,:)-out.centerz(1:p1)*out.r*breg;

    if rmsep==1
        rmse=rrmse(x,y,h,kmax,'RSIMPLS',0,k);  
        out.rmsep=rmse.rmsep;
    end
    
    % testing several output parameters
    if ~isnan(rmsep_value)
        out.rmsep = rmsep_value;
        options.rmsep = 1;
    end
    
    if ~isnan(rsquared_value)
        out.rcs = rsquared_value;
    end
    if ~isnan(rss_value)
        out.rcs = [out.rcs;sqrt(rss_value)];
    end
    if ~isnan(rmsecv_value)
        gammahalf = 0.5*sqrt(rss_value) + 0.5*rmsecv_value;
        out.rcs = [out.rcs;gammahalf;rmsecv_value];
        options.rmsecv = 1;
    end
    if any(isnan(rsquared_value)) && any(isnan(rss_value)) && any(isnan(rmsecv_value))
        out.rcs = 0;
    end

   %The output:
   out.T=out.t;
   out.weights.p=out.p;
   out.weights.r=out.r;
   out.kr=kr;
   out.h=h;
   out.alpha=alfa;
   out.slope=b;
   out.int=int;
   out.yhat=x*b+repmat(int,n,1);
   out.x=x;
   out.y=y;
   out.res=y-out.yhat;
   out.class='RSIMPLS';
   out.k=k;
   out.cov=robpcareg.cov;
    
   if ~st
       %calculation of robust distances
       %Score distance
       out.Tcov=robpcareg.sigma(1:k,1:k);
       out.Tcenter=robpcareg.center(1:k);
       out.sd=sqrt(mahalanobis(out.t,out.Tcenter,'cov',out.Tcov))';
       out.cutoff.sd=sqrt(chi2inv(0.975,k));
       %robust residual distance
       if q1==1
           out.resd=out.res/sqrt(out.cov);
       else
           out.resd=sqrt(mahalanobis(out.res,zeros(1,q1),'cov',out.cov))';
       end
       %robust orthogonal distances
       xtilde=out.t*out.p';
       Cdiff=out.xcentr-xtilde;
       for i=1:n
           out.od(i,1)=norm(Cdiff(i,:));
       end
       r=rank(x);
       if k~=r
           [m,s]=unimcd(out.od.^(2/3),out.h);
           out.cutoff.od = sqrt(norminv(0.975,m,s)^3);
       else
           out.cutoff.od=0;
       end
       out.cutoff.resd=sqrt(chi2inv(0.975,q1));

       %Computing flags
       out.flag.od=out.od<=out.cutoff.od;
       out.flag.resd=abs(out.resd)<=out.cutoff.resd;
       out.flag.all=(out.flag.od & out.flag.resd);


       %Multivariate Rsquared
       Yw=y(out.flag.all==1,:);
       cYw=mcenter(Yw);
       res=out.res(out.flag.all==1,:);
       out.rsquared=1-(det(res'*res)/det(cYw'*cYw));

       %Assigning output
       if options.rmsep~=1 && options.rmsecv~=1
           out.rmsep=0;
       end

       if classic
           resultclassic=csimpls(x,y,'k',k,'plots',0);
       else
           resultclassic=0;
       end
       result=struct('slope',{out.slope}, 'int',{out.int},'fitted',{out.yhat},'res',{out.res}, 'cov',{out.cov},...
           'T',{out.T}, 'weights', {out.weights},'Tcenter',{out.Tcenter},'Tcov',{out.Tcov},'rsquared',{out.rsquared},'rcs',{out.rcs},'rmsep',{out.rmsep},...
           'k',{out.k},'alpha',{out.alpha},'h',{out.h},'sd', {out.sd},'od',{out.od},...
           'resd',{out.resd},'cutoff',{out.cutoff},'flag',{out.flag},'class',{out.class},'classic',{resultclassic},'robpca',{out.robpca});
       if result.rcs==0
           result=rmfield(result,'rcs');
       end
       if result.rmsep==0
           result=rmfield(result,'rmsep');
       end
   else
       result=out;
   end
end

% Plots
try
    if plots && options.classic
        makeplot(result,'classic',1,'labsd',labsd,'labod',labod,'labresd',labresd)
    elseif plots
        makeplot(result,'labsd',labsd,'labod',labod,'labresd',labresd)
    end
catch %output must be given even if plots are interrupted
end
