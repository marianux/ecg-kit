function [P,T,L,r,centerX,cX]=classSVD(x)

%CLASSSVD performs the singular value decomposition of a matrix with more
% rows than columns (uses svd.m)
%
% Required input: 
%          x : data matrix of size n by p where n>p
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Sabine Verboven, Mia Hubert
% Last Update: 06/07/2004

[n,p]=size(x); 

if n==1
    error('The sample size is 1. No SVD can be performed.')
end
cX=mean(x);
centerX=x-repmat(cX,n,1); 
[U,S,loadings]=svd(centerX./sqrt(n-1),0); 
eigenvalues=diag(S).^2;
tol = max([n p])*eigenvalues(1)*eps;
r= sum(eigenvalues>tol);
L=eigenvalues(1:r);
P=loadings(:,1:r);
T=centerX*P; 

