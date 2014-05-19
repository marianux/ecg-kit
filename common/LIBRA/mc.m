function [result] = mc(x)

%MC calculates the medcouple, a robust measure of skewness.
%
% The medcouple is described in:
%    Brys, G., Hubert, M. and Struyf, A. (2004),
%    "A Robust Measure of Skewness",
%    Journal of Computational and Graphical Statistics,
%    13, 996-1017. 
%
% Required input arguments:
%    x : Data matrix.
%        Rows of x represent observations, and columns represent variables.  
%        Missing values (NaN's) and infinite values (Inf's) are allowed, since observations (rows) 
%        with missing or infinite values will automatically be excluded from the computations.
%   
% The output of MC is a vector containing the medcouple for each column
%     of the data matrix x.
%
% Example:
%    result = mc([chi2rnd(5,1000,1) trnd(3,1000,1)]);
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Guy Brys
% Last Update: 31/07/2007

if (nargin<1)
    error('No input arguments')
end
[n,p] = size(x);
if n==1
    x = x';
    p = 1;
end
x(sum(~isfinite(x),2)>0,:)=[];
n = size(x,1);
if (n>50000)
    error('When there are more than 50000 observations, the MC may be uncomputable due to memory limitations.')
elseif (n>100)
    result = mcc(x);
elseif n==0
    error('Due to missing values, the data matrix is empty and the MC can not be computed.')
else
    result = 0.5*(-mcc(repmat(max(x),size(x,1),1)-x)+mcc(x));
end
