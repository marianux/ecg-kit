function result = removeObsMcd(data,i,inputH0,inputFull,csteps);

%REMOVEOBSMCD is an auxiliary function to perform cross-validation with MCD 
% (see cvMcd.m).
% 
% The input:
%      data : the original data
%         i : the index of the observation that has to be removed.
%   inputH0 : a structure that contains the following fields:
%      inputH0.H0 : the optimal H subset based on the original data
%      inputH0.j : the index of the observation that is removed within H0.
%      inputH0.same
% inputFull : a structure that contains the following fields:
%      inputFull.T : the mean data(H0,:)
%      inputFull.S : the cov data(H0,:)
%    csteps : csteps.value = 1 (default) if csteps must performed on the updated cov and mu. else = 0.
%             csteps.number = the maximal number of csteps to perform. default = 20.        
%
% The output:
%   out.P_min_i      : the loadingvector after observation i is removed.
%   out.L_min_i      : the eigenvalues after observation i is removed.
%   out.mu_min_i     : the center of the data after observation i is removed.
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by S. Engelen

if nargin < 5
    csteps.value = 1;
    csteps.number = 2;
end

n = size(data,1);
p = size(data,2);
h = length(inputH0.H0);

Tfull  = inputFull.T;
Sfull = inputFull.S;

H0 = inputH0.H0;
j = inputH0.j;

data_min_i = removal(data,i,0);

[upS,upT] = updatecov(data(H0,:),Sfull,Tfull,0,[j],0);
mahdist = mahalanobis(data_min_i,upT,'cov',upS);
sortmahdist = sort(mahdist);
factor=sortmahdist(h-1)/chi2inv((h-1)/(n-1),p/2); % adapted factor

% performing c-steps on the updated mu and cov.
if csteps.value
    oldobj = det(upS);
    for noCsteps = 1:csteps.number
        [mahsort, indsort] = sort(mahdist);
        covcstep = cov(data_min_i(indsort(1:(h-1)),:));
        mucstep = mean(data_min_i(indsort(1:(h-1)),:));
        obj = det(covcstep);       
        if abs(oldobj - obj) > 1.e-12
            oldobj = obj;
        else
            break
        end
    end
end

TMCD = upT;
SMCD = factor*upS;
mahdist = mahdist/factor;
weights = (mahdist <= chi2inv(0.975,p));

% Reweighting:
[mu_min_i,S_min_i] = weightmecov(data_min_i,weights);
[PS1,LS1] = eig(S_min_i);
[P_min_i,L_min_i] = sortPL(PS1,LS1);

res.mu_min_i = mu_min_i;
res.L_min_i = L_min_i;
res.P_min_i = P_min_i;

result = res;

%---------------------------------------------------------------------------------------------
function [sortedP,sortedL] = sortPL(P,L)

% Sorts P and L for PCA, the columns are sorted by decreasing eigenvalues.

[sortedEw,indexsortedEw] = greatsort(diag(L));
sortedP = P(:,indexsortedEw);
sortedL = diag(sortedEw);