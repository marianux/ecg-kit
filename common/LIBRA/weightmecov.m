function [wmean,wcov]=weightmecov(data,weights)

%WEIGHTMECOV computes the reweighted mean and covariance matrix of multivariate data.
% 
% Required input arguments:
%      data : data matrix
%   weights : weights of the observations
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
%Written by Katrien Van Driessen and Carlos Lopez
%Implemented on 29 July 2005
%Last updated: 05 March 2007

n = size(data,1);
nvar = size(data,2);

if ~(isempty(find(weights<0)))
    error('The weights are negative');
end

if size(weights,1)==1
   weights=weights';
end

%Using sparse matrix; the expression is valid even for non-double data type
q=find(weights);
wmean=sum(spdiags(double(weights/sum(weights)),0,length(weights),length(weights))*double(data));
wcov=(data(q,:)-repmat(wmean,length(q),1))'*(data(q,:)-repmat(wmean,length(q),1))/(sum(weights.^2)-1);
