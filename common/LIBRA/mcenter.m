function [mcX,mX]=mcenter(x)

%MCENTER mean-centers the data matrix x columnwise
%
% Required input arguments:
%       x : Data matrix (n observations in rows, p variables in columns)
% 
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Sabine Verboven

[n,p]=size(x);
mX=mean(x);
mcX=x-repmat(mX,n,1);
