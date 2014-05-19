function result = cvRobpca(data,kmax,resrob,rawres,h,csteps)

%CVROBPCA calculates the robust cross-validated PRESS (predicted residual error sum of squares) curve
% for ROBPCA in a fast way. This curve can be used to make a selection of the optimal number of 
% components. The function is used in robpca.m. 
%
% Input arguments:
%    data : the whole data set
%    kmax : the maximal number of components to be considered.
%  resrob : result of robpca for k = kmax on the whole data set. 
%  rawres : Optional input parameter. If equal to 1 (default), then the R-PRESS is calculated
%           based on the orthogonal distances of the points, that are calculated without performing the 
%           MCD step within ROBPCA. 
%       h : the quantile used in ROBPCA.
%  csteps : optional : csteps.value : whether c-steps should be performed (default) or not.                        
%                      csteps.number: the number of c-steps that must be performed. 
%
% Output: 
%   result.press    : the R-PRESS values when the minimum is used to define the weights.
%   result.weights  : the minimum, median and smooth weights.
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Sanne Engelen 
% Last Update: 01/07/2004, 03/07/2006
% Last Revision: 03/07/2006


n = size(data,1);
p = size(data,2);
r = rank(data);
alfa = 0.75;
teller_if_lus = 0;

% some initialisations
Pk_min_i = [];
Lk_min_i = [];
muk_min_i = [];
Tik = [];

if nargin < 4
    rawres = 1;
    alfa=0.75;
    h = min(floor(2*floor((n+kmax+1)/2)-n+2*(n-floor((n+kmax+1)/2))*alfa),n);
    csteps.value = 1;
    csteps.number = 2;
elseif nargin == 4
    alfa=0.75;
    h = min(floor(2*floor((n+kmax+1)/2)-n+2*(n-floor((n+kmax+1)/2))*alfa),n);
    csteps.value = 1;
    csteps.number = 2;
elseif nargin == 5
    csteps.value = 1;
    csteps.number = 2;
else
    csteps.value = csteps.value;
    if csteps.value
        csteps.number = csteps.number;
    end
end

outWeights = weightscvRobpca(data,r,resrob,kmax,h);
result.weights = outWeights;
w_min = outWeights.w_min;

% inputH0
H0 = resrob.Hsubsets.H0;
same.value = 0;

for i = 1:n
    same.value = 0;
    if isempty(find(H0 == i))
        if teller_if_lus >= 1
            same.value = 1;
        end
        teller_if_lus = teller_if_lus + 1;
    end
    
    if ~rawres
        % inputH1
        H1 = resrob.Hsubsets.H1;
        
        % Hfreq 
        Hfreq = resrob.Hsubsets.Hfreq;
    
        Hsets = [H0;H1;Hfreq];
        factor = 1;
        Hsets_min_i = RemoveObsHsets(Hsets,i);
        res = removeObsRobpca(data,i,kmax,Hsets_min_i,same,factor,csteps);
    else
        data_min_i = removal(data,i,0);
        H0_min_i = RemoveObsHsets(H0,i);
        if ~same.value
            [Pk_min_i,TH0,LH0_min_i,rH0_min_i,centerX,muk_min_i]=kernelEVD(data_min_i(H0_min_i,:));
            res.Pk_min_i = Pk_min_i;
            res.muk_min_i = muk_min_i;
        else
            res = same.res;
            Pk_min_i = res.Pk_min_i;
            muk_min_i = res.muk_min_i;
        end
    end

    if isempty(find(H0 == i))
        same.res = res;
    end
    
    Pkmax_min_i = res.Pk_min_i;
    muk_min_i = res.muk_min_i;
    
    for k = 1:kmax
        xhoed_k(i,(k-1)*p + 1:k*p) = (data(i,:) -  muk_min_i)*Pkmax_min_i(:,1:k)*Pkmax_min_i(:,1:k)' + muk_min_i;
        if k~=r
            odk(i,k) = norm(data(i,:) - xhoed_k(i,(k-1)*p + 1:k*p));
        else
            odk(i,k) = 0;
        end
    end
end

for k = 1:kmax
    press_min(k) = 1/sum(w_min)*w_min*odk(:,k).^2;
end

result.press = press_min;
%-----------------------------------------------------------------------------------------------------------
function Hsets_min_i = RemoveObsHsets(Hsets,i)

% removes the right index from the $h$-subsets in Hsets to 
% obtain (h - 1)-subsets.
% every h-subset is put as a row in Hsets.
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
%-----------------------------------------------------------------------------------------------------------
function out = weightscvRobpca(data,r,resrob,kmax,h)

% computes the weights needed for the calculation of the R-PRESS.
%
% input:
%   data : the whole data set.
%      r : the rank of the data set.
%   resrob: the result of robpca on the whole data set for k = kmax.
%   kmax : the maximal number of components to be considered.
%   h : (n-h+1) measures the number of outliers the algorithm should 
%       resist. Any value between n/2 and n may be specified. 
%
% output: 
%   out.w_min    : the minimum weights


odk = [];
n = size(data,1);
p = size(data,2);

% Determine fixed weights:
for k = 1:kmax
    muk = resrob.M;
    xhoed_k(:,(k-1)*p + 1:k*p) = (data -  repmat(muk,n,1))*resrob.P(:,1:k)*resrob.P(:,1:k)' + repmat(muk,n,1);
    
    for i = 1:n
        % defining the od for each k
        if k~=r
            odk(i,k) = norm(data(i,:) - xhoed_k(i,(k-1)*p + 1:k*p));
        else
            odk(i,k) = 0;
        end
    end
    
    % defining weights for odk:
    if k~=r
        [m,s]=unimcd(odk(:,k).^(2/3),h);
        cutoff(k)=sqrt(norminv(0.975,m,s).^3);
        wod(:,k) = (odk(:,k) <= cutoff(k));
    else
        cutoff(k) = 0;
        wod(:,k) = 1;
    end
end

% determine the weights for every observation:
wk = wod;

if size(wk,1) == 1 || size(wk,2) == 1
    w_min = wk';
else
    w_min = min(wk,[],2)';
end

out.w_min = w_min;
out.odk = odk;
out.cutoff = cutoff;