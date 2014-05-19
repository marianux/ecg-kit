function [P,T,L,r,centerX,cX]=kernelEVD(X);

%KERNELEVD (kernel eigenvalue decomposition) performs a singular value decomposition
% of a matrix with more columns than rows.
%
% Required input: 
%      x : data matrix of size n by p where n < p (else classSVD is invoked)
%
% I/O: [P,T,L,r,centerX,cX]=kernelEVD(X);
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Sabine Verboven, Mia Hubert
% Last Update: 17/06/2003 

[n,p]=size(X);
if n > p 
    [P,T,L,r,centerX,cX]=classSVD(X);
else
    cX=mean(X);
    centerX=X-repmat(cX,n,1);
    [P,L]=eig(centerX*centerX'/(n-1));
    [L,I]=greatsort(diag(L));
    P=P(:,I);
    tol=n*max(L)*eps;
    r=sum(L>tol);
    L=L(1:r);
    loadings=(centerX/sqrt(n-1))'*P(:,1:r)*diag(1./sqrt(L)); 
    %normalizing loadings by dividing by the sqrt(eigenvalues)
    T=centerX*loadings;
    P=loadings;
end