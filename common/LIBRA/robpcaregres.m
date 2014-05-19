function robpcareg=robpcaregres(x,y,w,cutoffWeights)

%ROBPCAREGRES is a function for robust multivariate regression, based on
% the output from RSIMPLS and ROBPCA. 
%
% Required input arguments:
%            x  : Data matrix of the explanatory variables
%                 (n observations in rows, p variables in columns).
%            y  : Data matrix of the response variables
%                 (n observations in rows, q variables in columns).
%            w  : Weight vector, which is zero for the outliers, and one for the
%                 regular observations.
% cutoffWeights : the cutoff to decide on the weights. 
%
% I/O: robpcareg=robpcaregres(x,y,w,cutoffWeights);
% 
% The output is a structure containing
%
%   robpcareg.x       : The original explanatory variables
%   robpcareg.y       : The original respons(es)
%   robpcareg.coeffs  : The estimated regression coefficients
%   robpcareg.resids  : The residuals
%   robpcareg.cov     : Estimated variance-covariance matrix of the errors 
%   robpcareg.weights : Final weights
%   robpcareg.center  : Robust center of [x,y];
%   robpcareg.sigma   : Robust covariance estimate of [x,y]
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%       http://wis.kuleuven.be/stat/robust.html
%
% Written by Karlien Vanden Branden
% Last update: 05/12/2002. 

dat=[x,y];							
[n,m]=size(dat);
intercept=ones(n,1);  				
q=size(y,2);		
k=size(x,2);

if nargin < 4
    cutoffWeights = sqrt(chi2inv(0.975,q));
end

%calculation of reweighted covariance and center
rawclasscov=cov(dat(w==1,:)) ;			
rawclasscenter=mean(dat(w==1,:));  
rawcentx=rawclasscenter(1:k)';						
rawcenty=rawclasscenter((k+1):k+q)';					
rawsigmax=rawclasscov(1:k,1:k);
rawsigmay=rawclasscov(k+1:k+q,k+1:k+q);
rawsigmaxy=rawclasscov(1:k,(k+1):k+q);
rawsigmayx=rawsigmaxy';

%calculation of slope and intercept + residuals
rawbeta=[inv(rawsigmax)*rawsigmaxy; (rawcenty-(rawsigmayx*inv(rawsigmax)*rawcentx))'];
rawcovE=rawsigmay-rawbeta(1:k,1:q)'*rawsigmax*rawbeta(1:k,1:q);
rawfitted=[x,intercept]*rawbeta(1:(k+1),:);
rawresid=y-rawfitted;

%calculation of the reweighted weights based on the residuals of the reweighted beta-coefficients
rewweights=zeros(n,1);
for j=1:n
    if (sqrt(rawresid(j,1:q)*pinv(rawcovE)*rawresid(j,1:q)')) <= cutoffWeights
        rewweights(j)=1;
    end
end

%regression reweighting part
rewclasscov=cov(dat(rewweights==1,:)); 
rewclasscenter=mean(dat(rewweights==1,:));

rewcenterx=rewclasscenter(1:k)';
rewcentery=rewclasscenter((k+1):m)';

rewsigmax=rewclasscov(1:k,1:k);
rewsigmaxy=rewclasscov(1:k,(k+1):m);
rewsigmayx=rewsigmaxy';
rewsigmay=rewclasscov((k+1):m,(k+1):m);

rewbetarew=[inv(rewsigmax)*rewsigmaxy; (rewcentery-(rewsigmayx*inv(rewsigmax)*rewcenterx))'];
rewE2=rewsigmay-rewbetarew(1:k,1:q)'*rewsigmax*rewbetarew(1:k,1:q);
rewfittedrew=[x,intercept]*rewbetarew(1:(k+1),:);
rewresidrew=y-rewfittedrew;

%output
robpcareg.x=x;
robpcareg.y=y;
robpcareg.coeffs=rewbetarew;
robpcareg.cov=rewE2; %based on location/scale and regression reweighting
robpcareg.resids=rewresidrew;
robpcareg.weights=rewweights;
robpcareg.sigma=rawclasscov;
robpcareg.center=rawclasscenter;