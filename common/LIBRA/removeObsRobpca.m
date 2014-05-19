function res = removeObsRobpca(data,i,k,Hsets_min_i,same,factor_ind,csteps)

%REMOVEOBSROBPCA is an auxiliary function to perform cross-validation with ROBPCA, 
% RPCR, RSIMPLS, RSIMCA (see cvRobpca.m, cvRpcr.m, cvRsimpls.m).
%
% Input: 
%    data : the data set
%       i : the observation that is removed, index with respect to the whole data set.
%       k : the number of principal components that has to be calculated.
%  Hsets_min_i : contains H0_min_i, H1_min_i and Hfreq_min_i as first, second and third row respectively.
%                The h-subsets are implemented by means of indices of the observations in data_min_i, which is
%                is the original data set minus sample i. 
%  same : structure : 
%           same.value : indicates whether some part of this algorithm can be skipped (= 1) or not ( = 0). 
%           same.res : if same.value = 1, then some additional information is needed. 
%  factor : optional. default  = 0, then the original consistency factor is used.
%             Else factor = 1, the consistency factor is adapted to the kmax approach.
%  csteps : optional structure: 
%            csteps.value : 1 (default, then csteps are performed within robpca), 0 (no csteps are performed within robpca)
%            csteps.number : the number of csteps that need to be performed (default = 2).
%
% Output:
% res is the result structure. It contains:
%   Pk_min_i  : update of the loadingmatrix for a certain k when observation i is deleted.
%   Lk_min_i  : update of the eigenvalues for a certain k when observation i is deleted.
%   muk_min_i : update of the center for a certain k when observation i is deleted.
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at:
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by S. Engelen

rot = [];
center = [];
P1 = [];

n = size(data,1);
p = size(data,2);
h = size(Hsets_min_i,2) + 1;

if nargin < 6    
    factor_ind = 0;    
    csteps.number = 2;   
    csteps.value = 1;
elseif nargin == 6    
    csteps.number = 2;    
    csteps.value = 1;
end

data_min_i= removal(data,i,0);
mu0 = mean(data);
center = mu0;

H0_min_j = Hsets_min_i(1,:);
H1_min_t = Hsets_min_i(2,:);

if ~same.value   
    [PH0_min_i,TH0,LH0_min_i,rH0_min_i,centerX,mu1trafo]=kernelEVD(data_min_i(H0_min_j,:));  
    LH0_min_i = diag(LH0_min_i);   
    
    % calculation of the projection
    % res.T2tilde = T2tilde;   
    center = mu1trafo;    
    rot = PH0_min_i;    
    T2tilde = (data_min_i - repmat(mu1trafo,n-1,1))*PH0_min_i;  
    T2tilde = T2tilde(:,1:k);  
    rot = rot(:,1:k);       
    
    % defining the outputstructure res:  
    res.PH0_min_i = PH0_min_i; 
    res.LH0_min_i = LH0_min_i;   
    res.mu1trafo = mu1trafo;
else
    % defining the input structure input
    res = same.res;  
    PH0_min_i = res.PH0_min_i;  
    LH0_min_i = res.LH0_min_i;   
    mu1trafo = res.mu1trafo;    
    T2tilde = (data_min_i - repmat(mu1trafo,n-1,1))*PH0_min_i;  
    center = mu1trafo;    
    rot = PH0_min_i;  
    rot = rot(:,1:k);   
    T2tilde = T2tilde(:,1:k);
end

mah = mahalanobis(T2tilde,zeros(1,k),'invcov',1./diag(LH0_min_i(1:k,1:k)));
oldobj = prod(diag(LH0_min_i(1:k,1:k)));P4 = eye(k);

if csteps.value    
    oldobj = prod(diag(LH0_min_i(1:k,1:k)));  
    for j = 1:csteps.number     
        [mahsort, indsort] = sort(mah);      
        dataH1 = T2tilde(indsort(1:(h-1)),:);    
        [P,T,L,r3,Xm,clmX] = classSVD(dataH1);  
        obj = prod(L);      
        T2tilde = (T2tilde - repmat(clmX,n-1,1))*P;    
        center = center + clmX*rot';   
        rot = rot*P;      
        mah = mahalanobis(T2tilde,zeros(1,size(T2tilde,2)),'invcov',1./L);      
        P4 = P4*P;    
        if abs(oldobj - obj) > 1.e-12      
            oldobj = obj;    
        else
            break   
        end
    end
else
    obj = oldobj;
end

% extra reweighting:
if k~=r3  
    XRc= data_min_i-repmat(center,n-1,1);   
    Xtilde = T2tilde*rot'; 
    Rdiff = XRc-Xtilde;   
    for i=1:n     
        odh(i,1)=norm(Rdiff(i,:));   
    end
    [m,s]=unimcd(odh.^(2/3),h);    
    cutoffodh = sqrt(norminv(0.975,m,s).^3);     
    indexset = find(odh<=cutoffodh)';
    [P,Threw,Lrew,rrew,Xmrew,clmX]=kernelEVD(data_min_i(indexset,:));   
    center = clmX;   
    rot = P(:,1:k);
end
T2tilde = (data_min_i - repmat(center,n-1,1))*rot;

% Perform mcdcov on some H-subsets:
[res_min_i,raw_min_i] = mcdcov(T2tilde,'Hsets',Hsets_min_i,'h',h-1,'ntrial',250,'plots',0,'factor',1);
% perform the last part of ROBPCA:
if raw_min_i.objective < obj  
    z = res_min_i;
else
    sortmah = sort(mah);  
    if h==n      
        factor=1;  
    else
        if factor_ind == 1        
            factor = sortmah(h-1)/chi2inv((h-1)/size(T2tilde,1),size(T2tilde,2)/2); % adjusted factor.       
        else
            factor = sortmah(h-1)/chi2inv((h-1)/size(T2tilde,1),size(T2tilde,2)); % non adjusted factor    
        end
    end
    mah = mah/factor;  
    weights = mah <= chi2inv(0.975,size(T2tilde,2)); 
    [center_noMCD,cov_noMCD] = weightmecov(T2tilde,weights);  
    mah = mahalanobis(T2tilde,center_noMCD,'cov',cov_noMCD);  
    z.flag = (mah <= chi2inv(0.975,size(T2tilde,2)));   
    z.center = center_noMCD;  
    z.cov = cov_noMCD;
end

covf=z.cov;
centerf=z.center;

% The final PC:
[P6tilde,L6]=eig(covf);
[L6,I]=greatsort(diag(L6));
P6tilde=P6tilde(:,I);
T_min_i=(T2tilde-repmat(centerf,n-1,1))*P6tilde;

Pk_min_i=rot(:,1:k)*P6tilde;
Lk_min_i = L6;

res.Pk_min_i = Pk_min_i;
res.Lk_min_i = Lk_min_i;
res.muk_min_i = center + centerf*rot';
