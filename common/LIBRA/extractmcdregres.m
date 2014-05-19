function [Bk,intk,sigmayykmaxrew_k,sigmattkmaxrew_k,Hopt] = extractmcdregres(resMcd,T,y,kmax,n,q,k,h,cutoffWeights)

%EXTRACTMCDREGRES is an auxiliary function for cross-validation with RPCR and RSIMPLS. 
% It extracts mcd-regression coefficients for a certain k based on 
% the raw estimate of the slope and intercept of the regression with kmax coefficients. 

% Input arguments: 
%  resMCD : the results of the mcdregression on the data minus observation i for kmax components. 
%    kmax : the maximal number of components
%       T : the k scores from ROBPCA
%       k : the number of components. 
%
% Output : 
%    Bk : the slope for k components.
%  intk : the intercept for k components.
%  Hopt : remark that only for k < kmax a hopt is provided. For kmax Hopt needs to be calculated outside this function. 
%  sigmayykmaxrew_k and sigmattkmaxrew_k : estimates of covariance matrices that are needed in other functions. 
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Sanne Engelen 
% Last Update: 08/10/2004


geg = [T,y];
j = k;

mutkmaxraw_k = resMcd.Mu(1:j);
muykmaxraw_k = resMcd.Mu((kmax+1):(kmax + q));
sigmattkmaxraw_k = resMcd.Sigma(1:j,1:j);
sigmaytkmaxraw_k = resMcd.Sigma((kmax+1):(kmax + q),1:j);
sigmatykmaxraw_k = resMcd.Sigma(1:j,(kmax + 1):(kmax + q));
sigmayykmaxraw_k = resMcd.Sigma((kmax + 1):(kmax + q),(kmax + 1):(kmax + q));

% perform csteps:
Sigmaraw_k = [sigmattkmaxraw_k, sigmatykmaxraw_k ; sigmaytkmaxraw_k, sigmayykmaxraw_k];
Muraw_k = [mutkmaxraw_k;muykmaxraw_k]';

prevdet = 0;
if j ~= kmax
    geg_k = [T(:,1:j),y];
    for noCsteps = 1:10
        [sortMD,indsm] = sort(mahalanobis(geg_k,Muraw_k,'cov',Sigmaraw_k));
        Sigmaraw_k = cov(geg_k(indsm(1:h),:));
        Muraw_k = mean(geg_k(indsm(1:h),:));
        obj = det(Sigmaraw_k);
        if noCsteps >= 20 |abs(obj - prevdet) < 10^(-4)
            if noCsteps >= 20
                disp('no convergence in csteps')
            end
            break
        end
        prevdet = obj;
    end
    Hopt = indsm(1:h);
end

mukmax_k =  Muraw_k;
mutkmax_k = Muraw_k(1:j)';
muykmax_k = Muraw_k((j+1):(j + q))';
sigmakmax_k = Sigmaraw_k;
sigmattkmax_k = Sigmaraw_k(1:j,1:j);
sigmaytkmax_k = Sigmaraw_k((j+1):(j + q),1:j);
sigmatykmax_k = Sigmaraw_k(1:j,(j + 1):(j + q));
sigmayykmax_k = Sigmaraw_k((j + 1):(j + q),(j + 1):(j + q));

sigmattinvkmax_k = inv(sigmattkmax_k);
Braw_k = sigmattinvkmax_k*sigmatykmax_k;
intraw_k = (muykmax_k - (sigmaytkmax_k*sigmattinvkmax_k*mutkmax_k))';

% reweighting:
rewweights=zeros(n,1);
rewfitted=T(:,1:j)*Braw_k + repmat(intraw_k,n,1);
rewresid=y-rewfitted; 
rewE=sigmayykmax_k-Braw_k'*sigmattkmax_k*Braw_k;
for noObs=1:n
    if (sqrt(rewresid(noObs,1:q)*inv(rewE)*rewresid(noObs,1:q)')) <= cutoffWeights
        rewweights(noObs)=1;
    end
end

rewclasscov=cov(geg(rewweights==1,:));
rewclasscenter=mean(geg(rewweights==1,:));

mukmaxrew_k =  rewclasscenter;
mutkmaxrew_k = rewclasscenter(1:j)';
muykmaxrew_k = rewclasscenter((j+1):(j + q))';
sigmakmaxrew_k = rewclasscov;
sigmattkmaxrew_k = rewclasscov(1:j,1:j);
sigmaytkmaxrew_k = rewclasscov((j+1):(j + q),1:j);
sigmatykmaxrew_k = rewclasscov(1:j,(j + 1):(j + q));
sigmayykmaxrew_k = rewclasscov((j + 1):(j + q),(j + 1):(j + q));

sigmattinvkmaxrew_k = inv(sigmattkmaxrew_k);
Bk = sigmattinvkmaxrew_k*sigmatykmaxrew_k;
intk = (muykmaxrew_k - (sigmaytkmaxrew_k*sigmattinvkmaxrew_k*mutkmaxrew_k))';% 
rewclasscov=cov(geg(rewweights==1,:));
rewclasscenter=mean(geg(rewweights==1,:));

mukmaxrew_k =  rewclasscenter;
mutkmaxrew_k = rewclasscenter(1:j)';
muykmaxrew_k = rewclasscenter((j+1):(j + q))';
sigmakmaxrew_k = rewclasscov;
sigmattkmaxrew_k = rewclasscov(1:j,1:j);
sigmaytkmaxrew_k = rewclasscov((j+1):(j + q),1:j);
sigmatykmaxrew_k = rewclasscov(1:j,(j + 1):(j + q));
sigmayykmaxrew_k = rewclasscov((j + 1):(j + q),(j + 1):(j + q));

sigmattinvkmaxrew_k = inv(sigmattkmaxrew_k);
Bk = sigmattinvkmaxrew_k*sigmatykmaxrew_k;
intk = (muykmaxrew_k - (sigmaytkmaxrew_k*sigmattinvkmaxrew_k*mutkmaxrew_k))';