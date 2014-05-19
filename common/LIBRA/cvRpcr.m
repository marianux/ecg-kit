function result = cvRpcr(x,y,kmax,rmsecv,h,k)

%CVRPCR calculates the robust RMSECV (root mean squared error of cross-validation) curve
% for RPCR or the robust RMSEP (root mean squared error of prediction) value in a fast way. 
% The R-RMSECV curve can be used to make a selection of the optimal number of 
% components to include in the regression model. The function is used in rpcr.m. 
%
% Input arguments:
%          x      : the explanatory variables
%          y      : the response variables
%          kmax   : the maximal number of components to be considered in the model.
%          rmsecv : Optional. If equal to 1 (default), the rmsecv is computed. 
%                   Else, rmsecv = 0 and then the rss and R2 are computed. 
%          h      : the quantile used in RPCR
%          k      : optional, if equal to zero (default), the RMSECV is calculated. 
%                   If different from 0, the RMSEP value is calculated.
% Output:
%  If rmsecv = 1
%     result.rmsecv  : the r-rmsecv values (only if the RMSECV is computed)
%     result.rmsep   : the rmsep value (only if RMSEP is computed)
%     result.residu  : the residuals for every k = 1,...,kmax
%  result.outWeights : the fixed weights used in the robust version of the RMSE. 
%  result.R2         : the coefficient of determination for every k.
%  result.rss        : the RSS values for every k.
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Sanne Engelen 
% Last Update: 08/10/2004, 03/07/2006
% Last Revision: 03/07/2006

% some initialisations
n = size(x,1);
p = size(x,2);
q = size(y,2);
teller_if_lus = 0;
cutoffWeights = sqrt(chi2inv(0.975,q));

if nargin < 4
    alfa=0.75;
    h=floor(2*floor((n+kmax+2)/2)-n+2*(n-floor((n+kmax+2)/2))*alfa);
    k = 0;
    rmsecv = 1;
elseif nargin == 4
    alfa=0.75;
    h=floor(2*floor((n+kmax+2)/2)-n+2*(n-floor((n+kmax+2)/2))*alfa);
    k = 0;
elseif nargin == 5
    k = 0;
end

if q == 1   
    FixedWeights = weightscvLtsregres(x,y,kmax,h,k);
else
    FixedWeights = weightscvMcdregres(x,y,kmax,h,k,cutoffWeights);
end
kmax=FixedWeights.kmax;

if rmsecv 
    weights = FixedWeights.weightsk;
    resrob = FixedWeights.resrob;
    
    if k == 0
        w_min = FixedWeights.w_min;
    end
    if size(resrob.Hsubsets.H0,2)==1
        resrob.Hsubsets.H0=resrob.Hsubsets.H0';
    end
    Hsets = [FixedWeights.Hsets;resrob.Hsubsets.H0;resrob.Hsubsets.H1;resrob.Hsubsets.Hfreq];
    
    for i = 1:n
        disp(['observation ',num2str(i),' is left out'])
        x_min_i = removal(x,i,0);
        y_min_i = removal(y,i,0);
        
        same.value = 0;
        if isempty(find(resrob.Hsubsets.H0 == i))
            if teller_if_lus >= 1
                same.value = 1;
            end
            teller_if_lus = teller_if_lus + 1;
        end
        
        % constructing Hsets of right size: h - 1. 
        Hsets_min_i = RemoveObsHsets(Hsets,i);
        
        if k == 0
            res = removeObsRobpca(x,i,kmax,Hsets_min_i,same,1);
        else
            res = removeObsRobpca(x,i,k,Hsets_min_i,same);
        end
        
        if isempty(find(resrob.Hsubsets.H0 == i))
            same.res = res;
        end
        
        Prob_min_i = res.Pk_min_i;
        murob_min_i = res.muk_min_i;
        
        Tk_min_i = [];
        Lk_min_i = [];
        coeffk_min_i = [];
        
        if k == 0
            j_ind = kmax;
        else
            j_ind = k;
        end
        
        Tkmax_min_i = (x_min_i - repmat(murob_min_i,n-1,1))*Prob_min_i(:,1:j_ind);
        
        if q == 1
            [resLts2,rawLts2] = ltsregres(Tkmax_min_i,y_min_i,'plots',0,'Hsets',Hsets_min_i,'h',h-1);
        else
            if k == 0
                Mcdres = mcdcov([Tkmax_min_i,y_min_i],'h',h-1,'plots',0,'Hsets',Hsets_min_i);
                % the full mcdregres is not needed, we only need the results of mcdcov.
                resMcd2.Mu = Mcdres.center';
                resMcd2.Sigma = Mcdres.cov;
            end
        end
        
        if k == 0
            j_start1 = 1;
            j_end1 = kmax;
        else
            j_start1 = k;
            j_end1 = k;
        end
        
        for j = j_start1:j_end1
            if k == 0
                if q == 1
                    [Bk,intk,weightslts] = extractlts(rawLts2,x_min_i, y_min_i, murob_min_i, Prob_min_i,kmax,j,h-1,n-1);
                else
                    Tk_min_i = (x_min_i - repmat(murob_min_i,n-1,1))*Prob_min_i(:,1:j);
                    [Bk,intk,sigmayykmaxrew_k,sigmattkmaxrew_k] = extractmcdregres(resMcd2,Tk_min_i,y_min_i,kmax,n-1,q,j,h-1,cutoffWeights);
                    coeffk = [Bk;intk];
                end
            else
                if q == 1
                    Bk = resLts2.slope;
                    intk = resLts2.int;
                    weightslts = resLts2.flag;
                else
                    resMcdregres = mcdregres(Tkmax_min_i,y_min_i,'Hsets',Hsets_min_i,'plots',0,'h',h-1);
                    Bk = resMcdregres.slope;
                    intk = resMcdregres.int;
                    coeffk = [Bk;intk];
                end
            end
            finalB = Prob_min_i(:,1:j)*Bk;
            finalInt = intk - murob_min_i*finalB;
            yhat_min_i = x(i,:)*finalB + finalInt;
            residk_min_i(i,(j - 1)*q + 1:j*q) = y(i,:) - yhat_min_i;
            
            % calculation of the resd: 
            if q > 1
                if k == 0
                    rewE2=sigmayykmaxrew_k-coeffk(1:j,1:q)'*sigmattkmaxrew_k*coeffk(1:j,1:q);
                    cen=zeros(q,1);
                    resd(i,j)=sqrt(mahalanobis(residk_min_i(i,(j - 1)*q + 1:j*q),cen','cov',rewE2))'; %robust distances of residuals
                else
                    resd(i,j) = sqrt(mahalanobis(residk_min_i(i,(j - 1)*q + 1:j*q),zeros(q,1),'cov',resMcdregres.cov));
                end
            else
                if k == 0
                    scale=sqrt(sum(weightslts.*residk_min_i(i,(j - 1)*q + 1:j*q).^2)/(sum(weightslts)-1));
                    resd(i,j) = residk_min_i(i,(j-1)*q + 1:j*q)/scale;
                else
                    resd(i,j) = residk_min_i(i,(j-1)*q + 1:j*q)/resLts2.scale;
                end
            end
        end
    end
    
    if k == 0
        for j = 1:kmax
            resk = residk_min_i(:,(j-1)*q + 1:j*q);
            if q == 1
                rmsecv_min(j) = sqrt(1/sum(w_min)*w_min*(resk).^2);
            else
                rmsecv_min(j) = sqrt(1/sum(w_min)*w_min*(mean((resk').^2))');
            end
        end
        result.rmsecv = rmsecv_min;
        result.residu = residk_min_i;
    else
        weights = weights(:,k);
        if q == 1
            rmsep = sqrt(1/sum(weights)*weights'*(residk_min_i(:,k)).^2);
        else
            rmsep = sqrt(1/sum(weights)*weights'*(mean((residk_min_i(:,(k-1)*q + 1:k*q)').^2))');
        end
        result.rmsep = rmsep;
        result.residu = residk_min_i(:,(k-1)*q + 1:k*q);
    end
end

result.outWeights = FixedWeights;
result.rss =  FixedWeights.rss;
result.R2 = FixedWeights.R2;

%-------------------------------------------------------------------------
function Hsets_min_i = RemoveObsHsets(Hsets,i)

% Removes the correct index from the $h$-subsets in Hsets to 
% obtain (h - 1)-subsets.
% Every h-set is put as a row in Hsets.
% i is the index of the observation that is removed from the whole data.

for r = 1:size(Hsets,1)
    if ~isempty(find(Hsets(r,:)== i))
        Hsets_min_i(r,:) = removal(Hsets(r,:),0,find(Hsets(r,:) == i));
    else
        Hsets_min_i(r,:) = Hsets(r,1:(end-1));
    end
    
    for j = 1:length(Hsets_min_i(r,:))
        if Hsets_min_i(r,j) > i
            Hsets_min_i(r,j) = Hsets_min_i(r,j) - 1;
        end
    end
end
%------------------------------------------------------------------------
function out = weightscvLtsregres(x,y,kmax,h,k)

% Computes the weights needed in the R-RMSECV function.
%
% Input arguments: 
%   x    : the independent variables
%   y    : the response variables
%   kmax : the maximal number of components to be considered.
%   h    : the number of observations on which the calculations are based.
%   k    : if zero (default), the kmax approach is used (for RMSECV). 
%          Else robpca is calculated for a specific k (for RMSEP). 
%
% Output arguments:
%   if  k = 0 (RMSECV):
%       out.w_min    : the weights obtained by taking the minimum over all k 
%       out.weightsk : the weights for every observation and every k (n x kmax).
%       out.resrob   : the results of robpca on x for kmax components.
%   if k ~= 0 (RMSEP):
%       out.weightsk : the weights for a specific k
%       out.resrob   : the results of robpca on x for k components.
%   out.R2       : the weighted Rsquared for each value of k
%   out.rss      : the weighted rss for each value of k


n = size(x,1);
p = size(x,2);
q = size(y,2);

if nargin < 5
    k = 0;
end

if k == 0
    ResRobWhole = robpca(x,'plots',0,'kmax',kmax,'h',h,'k',kmax);
else
    ResRobWhole = robpca(x,'plots',0,'kmax',kmax,'h',h,'k',k);
end
kmax = ResRobWhole.kmax;
Trob = ResRobWhole.T;
Mrob = ResRobWhole.M;
Prob = ResRobWhole.P;
HsetsH0 = ResRobWhole.Hsubsets.H0;
HsetsH1 = ResRobWhole.Hsubsets.H1;
HsetsHfreqrob = ResRobWhole.Hsubsets.Hfreq;

[ResLtsWhole,RawLtsWhole]=ltsregres(Trob,y,'plots',0,'h',h);

HsetsHfreq = ResLtsWhole.Hsubsets.Hfreq;

if k == 0
    j_start = 1;
    j_end = kmax;
else
    j_start = k;
    j_end = k;
end

for j = j_start:j_end
    if k == 0
        if j ~= kmax
            [Bk,intk,weightslts,HsetsHopt(j,:)] = extractlts(RawLtsWhole,x,y,Mrob,Prob,kmax,j,h,n);
        else
            [Bk,intk,weightslts] = extractlts(RawLtsWhole,x,y,Mrob,Prob,kmax,j,h,n);
            HsetsHopt(kmax,:) =  ResLtsWhole.Hsubsets.Hopt';
        end
    else
        Bk = ResLtsWhole.slope;
        intk = ResLtsWhole.int;
        weightslts = ResLtsWhole.flag;
        HsetsHopt = ResLtsWhole.Hsubsets.Hopt';
    end
    coeffk=[Bk;intk];
    finalB = Prob(:,1:j)*Bk;
    finalInt = intk - Mrob*finalB;
    fitted=x*finalB + repmat(finalInt,n,1);
    residuals(:,j)=y-fitted;
    scale=sqrt(sum(weightslts.*residuals(:,j).^2)/(sum(weightslts)-1));
    s0=scale;
    weightsk(:,j)=abs(residuals(:,j)/scale)<=sqrt(chi2inv(0.975,1));
end
    
% Determining fixed weights.
if k == 0
    if kmax == 1
        w_min = weightsk';
    else
        w_min = min(weightsk');
    end
    out.w_min = w_min;

    yw = mean(y(w_min == 1,:));
    y2=(y-repmat(yw,n,1)).^2;
    R = residuals.^2;
    D=sum(y2(w_min==1,:));
    for j = 1:kmax
        R1=R(w_min==1,j); 
        rss(j) = sum(sum(R1));
        R2(j)=1-rss(j)/sum(D); 
    end
    out.rss = 1/sum(w_min)*rss;
    out.R2 = R2;
else
    yw = mean(y(weightsk(:,k) == 1,:));
    y2=(y-repmat(yw,n,1)).^2;
    R = residuals.^2;
    D=sum(y2(weightsk(:,k)==1,:));
    R1=R(weightsk(:,k)==1,:); 
    rss(k) = sum(sum(R1));
    R2(k)=1-rss(k)/sum(D); 
    out.rss = 1/sum(weightsk(:,k))*rss(k);
    out.R2 = R2(k);
end
out.weightsk = weightsk;
out.resrob = ResRobWhole;
if size(HsetsH0,2)==1
    HsetsH0=HsetsH0';
end
out.kmax=kmax;
out.Hsets = [HsetsHopt;HsetsHfreq;HsetsH0;HsetsH1;HsetsHfreqrob];
%-----------------------------------------------------------------------------
function [Bk,intk,weightslts,Hopt] = extractlts(rawltskmax, x_min_i, y_min_i, murob_min_i, Prob_min_i,kmax,k,h,n)

% This function extracts lts-regression coefficients for a certain k based on 
% the slope and intercept of the regression with kmax coefficients. 
%
% Input arguments: 
%  rawltskmax  : the raw results of ltsregression with kmax components. 
%  x_min_i     : the regressors without observation i.
%  y_min_i     : the dependent variable without observation i.
%  murob_min_i : the center estimated using robpca on the data minus observation i.
%  Prob_min_i  : the loadings matrix estimated using robpca on the data minus observation i. 
%  k           : the number of components
%
% Output arguments:
%  Bk         : the slope for k components. 
%  intk       : the intercept for k components.
%  weightslts : the weights needed for the reweighting. 
%  Hopt       : the optimal h-subset. Not availabe for kmax. 

j = k;
coeffkraw_min_i = rawltskmax.coefficients;
sloperaw = coeffkraw_min_i(1:j,:);
intraw = coeffkraw_min_i(end,:);
Tk_min_i = (x_min_i - repmat(murob_min_i,n,1))*Prob_min_i(:,1:j);

% perform csteps on the raw estimates of the parameters:
prevobj = 0;
if j ~= kmax
    for noCsteps = 1:10
        residu = y_min_i - Tk_min_i*sloperaw - repmat(intraw,n,1);
        [sortresid,indsr] = sort(residu.^2);
        obj = sum(sortresid(1:h));
        [Q,R] = qr([Tk_min_i(indsr(1:h),:) ones(h,1)],0);
        z = R\(Q'*y_min_i(indsr(1:h),1));
        sloperaw = z(1:j,:);
        intraw = z(j+1,:);
        if noCsteps >= 20 | abs(obj - prevobj) < 10^(-4)
           if noCsteps >= 20
               disp('no convergence in Csteps')
           end
           break
        end
        prevobj = obj;
    end
    Hopt = indsr(1:h)';
end

% reweighting:
factor =rawconsfactorlts(h,n);
residu = y_min_i - Tk_min_i*sloperaw - repmat(intraw,n,1);
[sortresid,indsr] = sort(residu.^2);
sh0 = sqrt(1/(h)*sum(sortresid(1:h)));
s0 = sh0*factor;
m = 2*(n)/asvarscalekwad((h),(n));
quantile = tinv(0.9875,m);
weightslts = abs(residu/s0) <= quantile;
[Q,R] = qr([Tk_min_i(weightslts == 1,:) ones(sum(weightslts),1)],0);
z = R\(Q'*y_min_i(weightslts == 1));
Bk = z(1:j,:);
intk = z(j+1,:);
%----------------------------------------------------------------------------------------------------------------------
function out = weightscvMcdregres(x,y,kmax,h,k,cutoffWeights)

% computes the weights needed in the R-PRESS function.
%
% Input arguments: 
%   x    : the independent variables
%   y    : the response variables
%   kmax : the maximal number of components to be considered.
%   h    : the number of observations on which the calculations are based.
%   k    : if zero (default), the kmax approach is used then (for RMSECV). 
%          Else robpca is calculated for a specific k (for RMSEP). 
%
% Output arguments:
%   if  k = 0 (RMSECV):
%       out.w_min    : the weights obtained by taking the minimum over all k 
%       out.weightsk : the weights for every observation and every k (n x kmax).
%       out.resrob   : the results of robpca on [X,Y] for kmax components.
%   if k ~= 0 (RMSEP):
%       out.weightsk : the weights for a specific k
%       out.resrob   : the results of robpca on [X,Y] for k components.
%   out.R2       : the weighted Rsquared for each value of k
%   out.rss      : the weighted rss for each value of k

n = size(x,1);
p = size(x,2);
q = size(y,2);

if nargin < 5
    k = 0;
end

if k == 0
    ResRobWhole = robpca(x,'plots',0,'k',kmax,'kmax',kmax,'h',h);
else
    ResRobWhole = robpca(x,'plots',0,'k',k,'kmax',kmax,'h',h);
end
kmax=ResRobWhole.kmax;
Trob = ResRobWhole.T;
Prob = ResRobWhole.P;
Mrob = ResRobWhole.M;
HsetsH0 = ResRobWhole.Hsubsets.H0;
HsetsH1 = ResRobWhole.Hsubsets.H1;
HsetsHfreqrob = ResRobWhole.Hsubsets.Hfreq;

ResMCDWhole = mcdregres(Trob,y,'h',h,'plots',0);
HsetsHfreq = ResMCDWhole.Hsubsets.Hfreq;
    
if k == 0
    for j = 1:kmax
        if j ~= kmax
            if k == 0
                Tk = (x - repmat(Mrob,n,1))*Prob(:,1:j);
                [Bk,intk,sigmayykmaxrew_k,sigmattkmaxrew_k,HsetsHopt(j,:)] = extractmcdregres(ResMCDWhole,Tk,y,kmax,n,q,j,h,cutoffWeights);
            else
                Tk = (x - repmat(Mrob,n,1))*Prob(:,1:j);
                [Bk,intk,sigmayykmaxrew_k,sigmattkmaxrew_k,HsetsHopt(j,:)] = extractmcdregres(ResMCDWhole,Tk,y,j,n,q,j,h,cutoffWeights);
            end
        else
            [Bk,intk,sigmayykmaxrew_k,sigmattkmaxrew_k] = extractmcdregres(ResMCDWhole,Trob,y,kmax,n,q,j,h,cutoffWeights);
            HsetsHopt(kmax,:) = ResMCDWhole.Hsubsets.Hopt;
        end
        coeffk = [Bk;intk];
        finalB = Prob(:,1:j)*Bk;
        finalInt = intk - Mrob*finalB;
        yhat = x*finalB + repmat(finalInt,n,1);
        residu(:,(j-1)*q+1:j*q) = y - yhat;
        
        cen=zeros(q,1)';
        cov=sigmayykmaxrew_k - coeffk(1:j,1:q)'*sigmattkmaxrew_k*coeffk(1:j,1:q);
        [nn,pp]=size(residu(:,(j-1)*q+1:j*q));
        resd = sqrt(mahalanobis(residu(:,(j-1)*q+1:j*q),cen,'cov',cov))';
        weightsk(:,j) = (abs(resd)<=cutoffWeights);
    end
else
    ResMCDWhole = mcdregres(Trob(:,1:k),y,'h',h,'plots',0);
    if size(ResMCDWhole.flag,2)~=1
        ResMCDWhole.flag = ResMCDWhole.flag';
    end
    
    weightsk(:,k) = ResMCDWhole.flag;
    
    HsetsHfreq = ResMCDWhole.Hsubsets.Hfreq;
    HsetsHopt = ResMCDWhole.Hsubsets.Hopt;
end


% Determining fixed weights.

if k == 0
    if kmax == 1
        w_min = weightsk';
    else
        w_min = min(weightsk');
    end
    out.w_min = w_min;
    
    yw = mean(y(w_min == 1,:));
    y2=(y-repmat(yw,n,1)).^2;
    R = residu.^2;
    D=sum(y2(w_min==1,:));
    for j = 1:kmax
        R1=R(w_min==1,(j-1)*q+1:j*q); 
        rss(j) = sum(sum(R1));
        R2(j)=1-rss(j)/sum(D); 
    end
    out.rss = 1/(q*sum(w_min))*rss;
    out.R2 = R2;
else
    yw = mean(y(weightsk(:,k) == 1,:));
    y2=(y-repmat(yw,n,1)).^2;
    R = ResMCDWhole.res.^2;
    D=sum(y2(weightsk(:,k)==1,:));
    R1=R(weightsk(:,k)==1,:); 
    rss = sum(sum(R1));
    R2=1-rss/sum(D); 
    out.rss = 1/(q*sum(weightsk(:,k)))*rss;
    out.R2 = R2;
end
out.weightsk = weightsk;
out.resrob = ResRobWhole;
if size(HsetsH0,2)==1
    HsetsH0=HsetsH0';
end
out.kmax=kmax;
out.Hsets = [HsetsHopt;HsetsHfreq;HsetsH0;HsetsH1;HsetsHfreqrob];

%---------------------------------------------------------------
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