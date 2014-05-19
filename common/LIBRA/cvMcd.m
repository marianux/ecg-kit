function result = cvMcd(data,kmax,resMCD,h)

%CVMCD calculates the robust cross-validated PRESS (predicted residual error sum of squares)
% curve for the MCD method in a fast way. 
%
% Input arguments: 
%   data   : the full data set
%   kmax   : the maximal number of components to be considered (mostly kmax = p).
%   resMCD : the result of mcdcov(data,'plots',0,'factor',1)
%   h      : the quantile used in MCD.
%
% output:
%   result.press     : vector of length kmax with the press values 
%   result.weights   : the weights for all observations
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Sanne Engelen 
% Last Update: 01/07/2004

% Some initialisations:
n = size(data,1);
p = size(data,2);
r = rank(data);
Pk = [];
Lk = [];
teller_if_lus = 0;

if nargin < 4
    alfa = 0.75;
    h=floor(2*floor((n+p+1)/2)-n+2*(n-floor((n+p+1)/2))*alfa);
end

outWeights = weightscvMcd(data,r,kmax,resMCD,h);

w_min = outWeights.w_min;

Hopt = resMCD.Hsubsets.Hopt;
inputH0.H0 = Hopt;
Tfull = mean(data(Hopt,:));
Sfull = cov(data(Hopt,:));

for i = 1:n 
    % deciding which index should be removed from H0.
    inputH0.same = 0;
    if isempty(find(inputH0.H0 == i))
        inputH0.j = h;
        if teller_if_lus >= 1
            inputH0.same = 1;
        end
        teller_if_lus = teller_if_lus + 1;
    else
        inputH0.j = find(inputH0.H0 == i);
    end
      
    % assigning the input variables:
    inputFull.T = Tfull;
    inputFull.S = Sfull;
    
    if ~inputH0.same
        res = removeObsMcd(data,i,inputH0,inputFull);
    end
    if (isempty(find(inputH0.H0 == i))) & (teller_if_lus == 1)
        resfixed = res;
    end
    if isempty(find(inputH0.H0 == i)) & (teller_if_lus ~= 1)
        res = resfixed;
    end

    P_min_i = res.P_min_i;
    L_min_i = res.L_min_i;
    mu_min_i = res.mu_min_i;
    
     for k = 1:kmax
        clear Pk Lk;        
        Pk = P_min_i(:,1:k);
        Lk = L_min_i(1:k,1:k);
        Xhoedk_min_i(i,(k-1)*p + 1:k*p) = (data(i,:) - mu_min_i)*Pk*Pk' + mu_min_i; 
        
        if k~=r
            odk(i,k) = norm(data(i,:) - Xhoedk_min_i(i,(k-1)*p + 1:k*p));
        else
            odk(i,k) = 0;
        end
    end
end

for k = 1:kmax
    press_min(k) = 1/sum(w_min)*w_min*odk(:,k).^2;
end

result.press = press_min;
result.weights = outWeights;

%----------------------------------------------------------------------------------
function out = weightscvMcd(data,r,kmax,resMCD,h)

% computes the weights used to calculate the robust PRESS values.
% 
% input: 
%   data   : the whole data
%       r  : the rank of the data
%   kmax   : the maximal number of components to be considered
%   resMCD : the result of mcdcov(data,'plots',0,'factor',1)
%   h      : the number of observations on which the computations are based.
%
% output: 
%   out.w_min    : the weights computed by taken the minimum over all k

% Some initialisations:
n = size(data,1);
p = size(data,2);
Pk = [];
Lk = [];
Tik = [];

[P,L] = eig(resMCD.cov);
[L,I] = greatsort(diag(L));
P = P(:,I);

for k = 1:kmax
    Pk = P(:,1:k);
    if h==n
        Lk=L(1:k);
    else
    Lk = chi2inv(h/n,k)/chi2inv(h/n,kmax/2)*L(1:k);% with correction for the factor
    end
    muk = resMCD.center;
    Xhoedk(:,(k-1)*p + 1:k*p) = (data - repmat(muk,n,1))*Pk*Pk' + repmat(muk,n,1); 
    Tk = (data - repmat(muk,n,1))*Pk;
    
    for i =1:n
        % defining the sd for the observation that is left:
        sdk(i,k) = sqrt(mahalanobis(Tk(i,:),zeros(1,k),'cov',diag(Lk)));
    
        % defining the od for the observation that is left:
        if k~=r
            odk(i,k) = norm(data(i,:) - Xhoedk(i,(k-1)*p + 1:k*p));
        else
            odk(i,k) = 1;
        end
    end
    
    % defining weights for odk and sdk:
    if k~=r
        [m,s]=unimcd(odk(:,k).^(2/3),h);
        cutoff(k)=sqrt(norminv(0.975,m,s).^3); 
        wod(:,k) = (odk(:,k) <= cutoff(k));
    else
        cutoff(k)= 0; 
        wod(:,k) = 1;
    end
    wsd(:,k) = (sdk(:,k) <= sqrt(chi2inv(0.975,k)));
end

% determine the weights for every observation:
wk = wsd & wod;

if size(wk,1) == 1 | size(wk,2) == 1
    w_min = wk';
else
    w_min = min(wk');
end

out.w_min = w_min;