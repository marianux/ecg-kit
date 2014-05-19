function result=predict(x,y,inmodel)

%PREDICT is a function which computes regression results for new data based on
% the output from a RPCR or RSIMPLS analysis.
%
% I/O: result=predict(x,y,inmodel)
%
% Required input arguments:
%   x: Data matrix of the explanatory variables of the new data
%   y: Data matrix of the response variables of the new data
%   inmodel: output from a RPCR or RSIMPLS model
%
% The output of PREDICT is a structure containing:
%
%   result.fitted  : Fitted response values of the new data
%   result.res     : Residuals of the new data 
%   result.sd      : Score distances of the new data 
%   result.od      : Orthogonal distances of the new data
%   result.resd    : Residual distances of the new data
%   result.cutoff  : Cutoff values for the score (result.cutoff.sd), orthogonal 
%                    (result.cutoff.od) and residual distances (result.cutoff.resd).
%                    We use 0.99 quantiles of the chi-squared distribution.            
%   result.flag    : The observations whose score distance is larger than 
%                    'result.cutoff.sd' receive a flag 'result.flag.sd' equal
%                    to zero (good leverage points). Otherwise 'result.flag.sd'
%                    is equal to one. 
%                    The components 'result.flag.od' and 'result.flag.resd' are
%                    defined analogously, and determine the orthogonal outliers, 
%                    resp. the bad leverage points/vertical outliers. 
%                    The observations with 'result.flag.od' and 'result.flag.resd'
%                    equal to zero, can be considered as calibration outliers and receive
%                    'result.flag.all' equal to zero. The regular observations and the good leverage
%                    points have 'result.flag.all' equal to one.
%   result.rmsep   : the root mean squared error of the non-outlying data 
%   result.class   : the name of the method used: 'RPCR' or 'RSIMPLS'
%   result.classic : if there was a classical output from rpcr or rsimpls,
%                    all computations based on these classical results are 
%                    done as well. Note that the root mean squared error is 
%                    then based on the non-outlying data from the robust
%                    analysis, in order to compare the rmsep values on the
%                    same data set.
%   
% Example:
%  result.pcr = rpcr(Xtrain,ytrain,'k',3,'plots',0);
%  result  = predict(Xtest,ytest,result.pcr);
%
% Written by S. Verboven, M. Hubert 
% Last Revision: 20/05/2008

X=x;
Y=y;
[n,p]=size(X);
[n,q]=size(Y);

%Fitted values and residuals of the new data
result.fitted=X*inmodel.slope + repmat(inmodel.int,n,1);
result.res=Y-result.fitted;

%Computing distances
if strcmp(inmodel.class,'RPCR')
    XRc=X-repmat(inmodel.robpca.M,n,1);
    result.T=XRc*inmodel.robpca.P;
    result.sd=sqrt(mahalanobis(result.T,zeros(size(result.T,2),1),'cov',inmodel.robpca.L))';
    Xtilde=result.T*inmodel.robpca.P';
    Rdiff=XRc-Xtilde;
    for i=1:n
        result.od(i,1)=norm(Rdiff(i,:));
    end

    if q >1
        cen=zeros(q,1)';
        result.resd=sqrt(mahalanobis(result.res,cen,'cov',inmodel.mcdreg.cov))';
    else
        result.resd=result.res/inmodel.lts.scale; 
    end
    %cutoffs
    quan=chi2inv(0.99,inmodel.k);
    result.cutoff.sd=quan^0.5;
    if inmodel.k~=rank(X)
        [m,s]=unimcd(inmodel.od.^(2/3),inmodel.h);
        result.cutoff.od = sqrt(norminv(0.99,m,s).^3);
    else
        result.cutoff.od=0;
    end
    result.cutoff.resd=sqrt(chi2inv(0.99,q));
elseif strcmp(inmodel.class,'RSIMPLS')
    XRc=X-repmat(inmodel.robpca.M(1:p),n,1);
    result.T=XRc*inmodel.weights.r;
    result.sd=sqrt(mahalanobis(result.T,inmodel.Tcenter,'cov',inmodel.Tcov))';
    Xtilde=result.T*inmodel.weights.p';
    Rdiff=XRc-Xtilde;
    for i=1:n
        result.od(i,1)=norm(Rdiff(i,:));
    end
    if q >1
        cen=zeros(q,1)';
        result.resd=sqrt(mahalanobis(result.res,cen,'cov',inmodel.cov))';
    else
        result.resd=result.res/sqrt(inmodel.cov);
    end
    %cutoffs
    quan=chi2inv(0.99,inmodel.k);
    result.cutoff.sd=quan^0.5;
    if inmodel.k~=rank(X)
        [m,s]=unimcd(inmodel.od.^(2/3),inmodel.h);
        result.cutoff.od = sqrt(norminv(0.99,m,s).^3);
    else
        result.cutoff.od=0;
    end
    result.cutoff.resd=sqrt(chi2inv(0.99,q));
end


% Defining flags 
result.flag.od=(result.od<=result.cutoff.od);
result.flag.sd=(result.sd<=result.cutoff.sd);
result.flag.resd=(abs(result.resd)<=result.cutoff.resd);
result.flag.all=result.flag.od & result.flag.resd;

% Rmsep based on data with 'result.flag.resd = 1'
N=sum(result.flag.resd);
if q>1
    result.rmsep=sqrt(1/(N*q)*sum(sum(result.res(result.flag.resd==1).^2,2)));
else
    result.rmsep=sqrt(1/N*sum((result.res(result.flag.resd==1)).^2));
end
result.class=inmodel.class;

%In case the classical output is also given
if isstruct(inmodel.classic) && strcmp(inmodel.class,'RPCR')
    % fitted values, residuals, distances
    result.classic.fitted=X*inmodel.classic.slope + repmat(inmodel.classic.int,n,1);
    result.classic.res=Y-result.classic.fitted; 
    XRc=X-repmat(inmodel.classic.cpca.M,n,1);
    result.classic.T=XRc*inmodel.classic.cpca.P;
    result.classic.sd=sqrt(mahalanobis(result.classic.T,zeros(size(result.classic.T,2),1),'cov',inmodel.classic.cpca.L))';
    Xtilde=result.classic.T*inmodel.classic.cpca.P';
    Rdiff=XRc-Xtilde;
    for i=1:n
        result.classic.od(i,1)=norm(Rdiff(i,:));
    end
    if q >1
        cen=zeros(q,1)';
        result.classic.resd=sqrt(mahalanobis(result.classic.res,cen,'cov',inmodel.classic.cov))';
    else
        result.classic.resd=result.classic.res/sqrt(inmodel.classic.cov);
    end
    % cutoffs and flags
    quan=chi2inv(0.99,inmodel.k);
    result.classic.cutoff.sd=quan^0.5;
    if inmodel.k~=rank(X)
        [m,s]=unimcd(inmodel.classic.od.^(2/3),inmodel.h);
        result.classic.cutoff.od = sqrt(norminv(0.99,m,s).^3);
    else
        result.classic.cutoff.od=0;
    end
    result.classic.cutoff.resd=sqrt(chi2inv(0.99,q));
    result.classic.flag.od=(result.classic.od<=result.classic.cutoff.od);
    result.classic.flag.sd=(result.classic.sd<=result.classic.cutoff.sd);
    result.classic.flag.resd=(abs(result.classic.resd)<=result.classic.cutoff.resd);
    result.classic.flag.all=result.classic.flag.od & result.classic.flag.resd;
    % classical rmsep of the good observations (from the robust analysis)
    N=sum(result.flag.resd);
    if q>1
        result.classic.rmsep=sqrt(1/(N*q)*sum(sum(result.classic.res(result.flag.resd==1).^2,2)));
    else
        result.classic.rmsep=sqrt(1/N*sum((result.classic.res(result.flag.resd==1)).^2));
    end
    result.classic.class=inmodel.classic.class;
elseif isstruct(inmodel.classic) && strcmp(inmodel.class,'RSIMPLS')
    % fitted values, residuals, distances
    result.classic.fitted=X*inmodel.classic.slope + repmat(inmodel.classic.int,n,1);
    result.classic.res=Y-result.classic.fitted;
    XRc=X-repmat(inmodel.classic.M(1:p),size(X,1),1);
    result.classic.T=XRc*inmodel.classic.weights.r;
    result.classic.sd=sqrt(mahalanobis(result.classic.T,zeros(size(result.classic.T,2),1),'cov',inmodel.classic.Tcov))';
    Xtilde=result.classic.T*inmodel.classic.weights.p';
    Rdiff=XRc-Xtilde;
    for i=1:n
        result.classic.od(i,1)=norm(Rdiff(i,:));
    end
    if q >1
        cen=zeros(q,1)';
        result.classic.resd=sqrt(mahalanobis(result.classic.res,cen,'cov',inmodel.classic.cov))';
    else
        result.classic.resd=result.classic.res/sqrt(inmodel.classic.cov);
    end
    % cutoffs and flags
    quan=chi2inv(0.99,inmodel.k);
    result.classic.cutoff.sd=quan^0.5;
    if inmodel.k~=rank(X)
        [m,s]=unimcd(inmodel.classic.od.^(2/3),inmodel.h);
        result.classic.cutoff.od = sqrt(norminv(0.99,m,s).^3);
    else
        result.classic.cutoff.od=0;
    end
    result.classic.cutoff.resd=sqrt(chi2inv(0.99,q));
    result.classic.flag.od=(result.classic.od<=result.classic.cutoff.od);
    result.classic.flag.sd=(result.classic.sd<=result.classic.cutoff.sd);
    result.classic.flag.resd=(abs(result.classic.resd)<=result.classic.cutoff.resd);
    result.classic.flag.all=result.classic.flag.od & result.classic.flag.resd;
    % classical rmsep of the good observations (from the robust analysis)
    N=sum(result.flag.resd);
    if q>1
        result.classic.rmsep=sqrt(1/(N*q)*sum(sum(result.classic.res(result.flag.resd==1).^2,2)));
    else
        result.classic.rmsep=sqrt(1/N*sum((result.classic.res(result.flag.resd==1)).^2));
    end
    result.classic.class=inmodel.classic.class; 
end


