function result = cvRsimpls(x,y,kmax,rmsecv,h,k)

%CVRIMPLS calculates the robust RMSECV (root mean squared error of cross-validation) curve
% for RSIMPLS or the robust RMSEP(root mean squared error of prediction) value in a fast way. 
% The R-RMSECV curve can be used to make a selection of the optimal number of 
% components to include in the regression model. The function is used in rsimpls.m. 
%
% Input arguments:
%    x          : the explanatory variables
%    y          : the rsponse variables
%    kmax       : maximal number of variables to include in the model.
%    rmsecv     : Optional. If equal to 1 (default), the rmsecv is computed. 
%                 Else, rmsecv = 0 and then the rss and R2 are computed. 
%    h          : optional input argument.
%    k          : optional input argument. If k = 0 (default) then RMSECV is calculated. Else the RMSEP wil be computed.
%
% Output: 
%  if RMSECV is computed:
%   result.rmsecv     : the R-RMSECV values (obtained with the minimum weights).
%  if RMSEP is computed: 
%   result.rmsep      : the R-RMSEP values
%   result.rss        : the RSS values for every k.
%   result.R2         : the coefficient of determination for every k.
%   result.residu     : the residuals for every k = 1,...,kmax
%   result.outWeights : the weights used to compute the robust R-RMSEP values.
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
r = rank(x);
rz = rank([x,y]);
teller_if_lus = 0;
cutoffWeights = sqrt(chi2inv(0.975,q));

if nargin < 4
    alfa = 0.75;
    kmaxr=min([kmax+q,rz]);
    h=floor(2*floor((n+kmaxr+1)/2)-n+2*(n-floor((n+kmaxr+1)/2))*alfa);
    k = 0;
    rmsecv = 1;
elseif nargin == 4
    alfa = 0.75;
    kmaxr=min([kmax+q,rz]);
    h=floor(2*floor((n+kmaxr+1)/2)-n+2*(n-floor((n+kmaxr+1)/2))*alfa);
    k = 0;
elseif nargin == 5
    k = 0;
end

outWeights = weightcvRsimpls(x,y,kmax,h,k,cutoffWeights);

if rmsecv
    if k == 0
        w_min = outWeights.w_min;
    end
    rob = outWeights.ResRob;
    
    % Assigning the input variables
    if size(rob.Hsubsets.H0,2)==1
        rob.Hsubsets.H0=rob.Hsubsets.H0';
    end
    Hsets = [rob.Hsubsets.H0;rob.Hsubsets.H1;rob.Hsubsets.Hfreq];
    same.value = 0;
    data = [x,y];
    
    for i = 1:n
        disp(['observation ',num2str(i),' is left out'])
        X_min_i = removal(x,i,0);
        Y_min_i = removal(y,i,0);
        data_min_i = removal(data,i,0);
        
        same.value = 0;
        if isempty(find(rob.Hsubsets.H0 == i))
            if teller_if_lus >= 1
                same.value = 1;
            end
            teller_if_lus = teller_if_lus + 1;
        end
        
        % constructing Hsets of right size: h - 1. 
        Hsets_min_i = RemoveObsHsets(Hsets,i);
        
        if k == 0
            res = removeObsRobpca(data,i,kmax + q,Hsets_min_i,same,1);
        else
            res = removeObsRobpca(data,i,k + q,Hsets_min_i,same);
        end
        
        if isempty(find(rob.Hsubsets.H0 == i))
            same.res = res;
        end
        
        Prob_min_i = res.Pk_min_i;
        Lrob_min_i = res.Lk_min_i;
        murob_min_i = res.muk_min_i;
        Trob_min_i = (data_min_i - repmat(murob_min_i,n-1,1))*Prob_min_i;
        
        % Computing weights corresponding with the ROBPCA results. 
        sdrob_min_i = sqrt(mahalanobis(Trob_min_i,zeros(1,size(Trob_min_i,2)),'invcov',1./Lrob_min_i))';
        
        if k == 0
            cutoff.sd=sqrt(chi2inv(0.975,kmax));
        else
            cutoff.sd = sqrt(chi2inv(0.975,k));
        end
        
        % Orthogonal distances to robust PCA subspace
        XRc=data_min_i-repmat(murob_min_i,n-1,1);
        Xtilde=Trob_min_i*Prob_min_i';
        Rdiff=XRc-Xtilde;
        for j=1:(n-1)
            odrob_min_i(j,:)=norm(Rdiff(j,:));
        end
        % Robust cutoff-value for the orthogonal distance
        if k == 0
            test_k = kmax;
        else
            test_k = k;
        end
        
        if test_k~=r
            [m,s]=unimcd(odrob_min_i.^(2/3),h);
            cutoff.od = sqrt(norminv(0.975,m,s).^3); 
            wrob_min_i = (odrob_min_i<=cutoff.od)&(sdrob_min_i<=cutoff.sd);
        else
            cutoff.od=0;
            wrob_min_i = (sdrob_min_i<=cutoff.sd);
        end
        
        % start the deflation:
        xycentr_min_i = [];
        sigmax_min_i = [];
        xcentr_min_i = [];
        sigmaxy_min_i = [];
        sigmayx_min_i = [];
        
        xycentr_min_i = murob_min_i;
        sigmax_min_i = Prob_min_i(1:p,:)*diag(Lrob_min_i)*Prob_min_i(1:p,:)';
        xcentr_min_i = X_min_i - repmat(murob_min_i(1:p),n-1,1);
        sigmaxy_min_i = Prob_min_i(1:p,:)*diag(Lrob_min_i)*Prob_min_i(p+1:p+q,:)';
        sigmayx_min_i = sigmaxy_min_i';
        
        % calculation of the scores.
        nScores = 1;
        R_min_i = [];
        T_min_i = [];
        P_min_i = [];
        V_min_i = [];
        
        if k == 0
            countScores = kmax;
        else
            countScores = k;
        end
        
        while nScores <= countScores
            if q == 1
                qq_min_i = 1;
            else
                [QQ,LL] = eig(sigmayx_min_i*sigmaxy_min_i);
                [LL,I] = greatsort(diag(LL));
                qq_min_i = QQ(:,I(1));
            end
            
            rr_min_i = sigmaxy_min_i*qq_min_i;
            rr_min_i = rr_min_i/norm(rr_min_i);
            tt_min_i = xcentr_min_i*rr_min_i;
            pp_min_i = sigmax_min_i*rr_min_i/(rr_min_i'*sigmax_min_i*rr_min_i);
            vv_min_i = pp_min_i;
            
            if nScores > 1
                vv_min_i = vv_min_i - V_min_i*(V_min_i'*pp_min_i);
            end
            vv_min_i = vv_min_i./norm(vv_min_i);
            sigmaxy_min_i = sigmaxy_min_i - vv_min_i*(vv_min_i'*sigmaxy_min_i);
            
            V_min_i(:,nScores) = vv_min_i;
            T_min_i(:,nScores) = tt_min_i;
            R_min_i(:,nScores) = rr_min_i;
            P_min_i(:,nScores) = pp_min_i;
            
            nScores = nScores + 1;
        end
        
        if k == 0
            outRegr = runRegr(x,y,i,T_min_i,Y_min_i,R_min_i,murob_min_i,wrob_min_i,kmax,cutoffWeights);
            for j = 1:kmax
                Tk_min_i = T_min_i(:,1:j);
                geg = [Tk_min_i,Y_min_i];
                
                outRegr.Mu = outRegr.center';
                outRegr.Sigma = outRegr.sigma;
                [Bk,intk,sigmayykmaxrew_k,sigmattkmaxrew_k] = extractmcdregres(outRegr,Tk_min_i,Y_min_i,kmax,n-1,q,j,h-1,cutoffWeights);
                coeffk = [Bk;intk];
                b_min_i = R_min_i(:,1:j)*coeffk(1:j,:);
                int_min_i = coeffk(j+1,:)  - murob_min_i(1:p)*R_min_i(:,1:j)*coeffk(1:j,:);
                Yhat_min_i = x(i,:)*b_min_i + int_min_i;
                resid_min_i(i,(j-1)*q + 1:j*q) = y(i,:) - Yhat_min_i;
                
                % calculation of the resd: 
                rewE2=sigmayykmaxrew_k- coeffk(1:j,1:q)'*sigmattkmaxrew_k*coeffk(1:j,1:q);
                
                if q > 1
                    cov = rewE2;
                    cen=zeros(q,1);
                    resd(i,j)=sqrt(mahalanobis(resid_min_i(i,(j-1)*q + 1:j*q),cen','cov',cov))'; %robust distances of residuals
                else
                    scale = sqrt(rewE2);
                    resd(i,j) = resid_min_i(i,(j-1)*q + 1:j*q)/scale;
                end
            end
        else
            outRegr = runRegr(x,y,i,T_min_i,Y_min_i,R_min_i,murob_min_i,wrob_min_i,k,cutoffWeights);   
            resid_min_i(i,:) = outRegr.resid_min_i;
            resd(i,:) = outRegr.resd';
        end
    end
    
    if k == 0
        for j = 1:kmax
            resk = resid_min_i(:,(j-1)*q + 1:j*q);
            if q == 1
                rmsecv(j) = sqrt(1/sum(w_min)*w_min*(resk).^2);
            else
                rmsecv(j) = sqrt(1/sum(w_min)*w_min*(mean((resk').^2))');
            end
        end
        result.rmsecv = rmsecv;
        result.residu = resid_min_i;
    else
        weights = outWeights.weightsk;
        if q == 1
            rmsep = sqrt(1/sum(weights)*weights'*(resid_min_i).^2);
        else
            rmsep = sqrt(1/sum(weights)*weights'*(mean((resid_min_i').^2))');
        end
        result.rmsep = rmsep;
        result.residu = resid_min_i;
    end
end

result.outWeights = outWeights;
result.rss = outWeights.rss;
result.R2 = outWeights.R2;

%------------------------------------------------------------------
function out = runRegr(x,y,i,T_min_i,Y_min_i,R_min_i,mukmax_min_i,wkmax_min_i,k,cutoffWeights)

[n,p] = size(x);
[n,q] = size(y);

% perform the robpca regression:
breg = [];
b_min_i = [];
int_min_i= [];
Yhat_min_i = [];
robpcareg = robpcaregres(T_min_i(:,1:k),Y_min_i,wkmax_min_i',cutoffWeights);
out.center = robpcareg.center;
out.sigma = robpcareg.sigma;
breg = robpcareg.coeffs(1:k,:);
b_min_i = R_min_i(:,1:k)*breg;
int_min_i = robpcareg.coeffs(k+1,:) - mukmax_min_i(1:p)*R_min_i(:,1:k)*breg;
Yhat_min_i = x(i,:)*b_min_i + int_min_i;
resid_min_i = y(i,:) - Yhat_min_i;

% calculation of the resd: 
if q > 1
    cov = robpcareg.cov;
    cen=zeros(q,1);
    resd=sqrt(mahalanobis(resid_min_i,cen','cov',cov))'; %robust distances of residuals
else
    scale = sqrt(robpcareg.cov);
    resd = resid_min_i/scale;
end

out.resid_min_i = resid_min_i;
out.resd = resd;
%----------------------------------------------------------------------------------------
function out = weightcvRsimpls(x,y,kmax,h,k,cutoffWeights)

% Computes the weights for the robust RMSECV/RMSEP values.
%
% input:
%   x    : the independent variables.
%   y    : the response variables.
%   kmax : the maximal number of components to be considered.
%   h    : the number of observations on which the calculations are based.
%   k    : if equal to zero, robpca is performed on kmax components (case RMSECV). (default). 
%          Else, robpca is performed for a certain number of components (case RMSEP)
%
% output: 
%   out.w_min    : the weights obtained by taking the minimum over all k
%   out.weightsk : the weights for all observations and all k (n x kmax)
%   out.resrob   : the results of robpca on [x,y].
%   out.R2       : the weighted Rsquared for each value of k
%   out.rss      : the weighted rss for each value of k

n = size(x,1);
p = size(x,2);
q = size(y,2);
r = rank(x);

if nargin < 5
    k = 0;
end

if k == 0
    ResRob = robpca([x,y],'plots',0,'k',kmax + q,'kmax',kmax + q,'h',h);
else
    ResRob = robpca([x,y],'plots',0,'k',k + q,'kmax',kmax + q,'h',h);
end

Trob = ResRob.T;
Prob = ResRob.P;
Lrob = ResRob.L;
murob = ResRob.M;
wrob = ResRob.flag.all;

%deflation

xycentr = [];
sigmax = [];
xcentr = [];
sigmaxy = [];
sigmayx = [];

xycentr = murob;
sigmax = Prob(1:p,:)*diag(Lrob)*Prob(1:p,:)';
xcentr = x - repmat(murob(1:p),n,1);
sigmaxy = Prob(1:p,:)*diag(Lrob)*Prob(p+1:p+q,:)';
sigmayx = sigmaxy';

% calculation of the scores.
nScores = 1;
R = [];
T = [];
P = [];
V = [];

if k == 0
    countScores = kmax;
else
    countScores = k;
end

while nScores <= countScores
    if q == 1
        qq = 1;
    else
        [QQ,LL] = eig(sigmayx*sigmaxy);
        [LL,I] = greatsort(diag(LL));
        qq = QQ(:,I(1));
    end
    
    rr = sigmaxy*qq;
    rr = rr/norm(rr);
    tt = xcentr*rr;
    pp = sigmax*rr/(rr'*sigmax*rr);
    vv = pp;
    
    if nScores > 1
        vv = vv - V*(V'*pp);
    end
    vv = vv./norm(vv);
    sigmaxy = sigmaxy - vv*(vv'*sigmaxy);
    
    V(:,nScores) = vv;
    T(:,nScores) = tt;
    R(:,nScores) = rr;
    P(:,nScores) = pp;
    
    nScores = nScores + 1;
end

if k == 0
    outRobRegr = robpcaregres(T(:,1:kmax),y,wrob');
    
    for j = 1:kmax
        outRobRegr.Mu = outRobRegr.center';
        outRobRegr.Sigma = outRobRegr.sigma;
        [Bk,intk,sigmayykmaxrew_k,sigmattkmaxrew_k] = extractmcdregres(outRobRegr,T(:,1:j),y,kmax,n,q,j,h,cutoffWeights);
        coeffk = [Bk;intk];
        finalB = R(:,1:j)*coeffk(1:j,:);
        finalInt = coeffk(j+1,:)  - murob(1:p)*R(:,1:j)*coeffk(1:j,:);
        Yhat = x*finalB + repmat(finalInt,n,1);
        resid(:,(j-1)*q+1:j*q) = y - Yhat;
        
        % calculation of the rd: 
        rewE2=sigmayykmaxrew_k- coeffk(1:j,1:q)'*sigmattkmaxrew_k*coeffk(1:j,1:q);
        
        if q > 1
            cov = rewE2;
            cen=zeros(q,1);
            resd = sqrt(mahalanobis(resid(:,(j-1)*q+1:j*q),cen','cov',cov))'; %robust distances of residuals
            weightsk(:,j) = (abs(resd)<=cutoffWeights);
        else
            scale = sqrt(rewE2);
            resd = resid(:,(j-1)*q+1:j*q)/scale;
            weightsk(:,j) = (abs(resd)<=cutoffWeights);
        end
    end
else
    outRobRegr = robpcaregres(T(:,1:k),y,wrob');
    
    % robust residual distance:
    if q==1
        resd=outRobRegr.resids/sqrt(outRobRegr.cov); 
    else
        resd=sqrt(mahalanobis(outRobRegr.resids,zeros(1,q),'cov',outRobRegr.cov))';
    end
    
    weightsk = (abs(resd)<=cutoffWeights);
end

if k == 0
    if kmax == 1
        w_min = weightsk';
    else
        w_min = min(weightsk');
    end
    
    out.w_min = w_min;
    out.weightsk = weightsk;
    
    yw = mean(y(w_min == 1,:));
    y2=(y-repmat(yw,n,1)).^2;
    R = resid.^2;
    D=sum(y2(w_min==1,:));
    for j = 1:kmax
        R1=R(w_min==1,(j-1)*q+1:j*q); 
        rss(j) = sum(sum(R1));
        R2(j)=1-rss(j)/sum(D); 
    end
    
    out.rss = 1/(q*sum(w_min))*rss;
    out.R2 = R2;
else
    out.weightsk = weightsk;
    
    s=sum(weightsk);
    yw=sum(y(weightsk==1,:))/s;  
    y2=(y-repmat(yw,n,1)).^2;
    R = outRobRegr.resids.^2;
    D=sum(y2(weightsk==1,:));
    R1=R(weightsk==1,:); 
    rss = sum(sum(R1));
    R2=1-rss/sum(D); 
    
    out.rss = 1/(q*sum(weightsk))*rss;
    out.R2 = R2;
end
out.ResRob = ResRob;
%---------------------------------------------------------------------------------------
function Hsets_min_i = RemoveObsHsets(Hsets,i)

% removes the right index from the $h$-subsets in Hsets to 
% obtain (h - 1)-subsets.
% every h-set is put as a row in Hsets.
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