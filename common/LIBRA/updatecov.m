function [upcovar,upmu] = updatecov(data,covar,mu,addrow,delrow,option)

% This function updates the mean and covariance matrix of the full data 
% after adding or removing observations.
%
% Input arguments: 
%              data : the original data
%             covar : the covariance matrix of the data
%                mu : the center of the data
%            addrow : matrix with contains the measurements of the data points that are added.
%            delrow : indices of the data points that are removed from the data set.
%            option : indicates whether you want to remove or add an observation.
%                     0 : remove
%                     1 : add
%                     2 : change a row with another vector. The row that had to be changed has to be put in
%                        delrow. The vector that needs to be add has to be put in addrow.

n = size(data,1);

if option == 0
    for i = 1:length(delrow)
       j = delrow(i);
       mu = (n-i+1)/(n - i)*(mu - 1/(n - i + 1)*data(j,:));
       covar = (n - i)/(n - i -1)*covar - (n - i)/((n - i + 1)*(n-i-1))*((mu - data(j,:))'*(mu - data(j,:)));
    end
    upmu = mu; 
    upcovar = covar;
elseif option == 1
    for i = 1:size(addrow,1)
       covar = (n+i-2)/(n +i-1)*covar + 1/(n+i)*((mu - addrow(i,:))'*(mu - addrow(i,:)));
       mu = (n+i-1)/(n+i)*mu + 1/(n+i)*addrow(i,:);
    end
    upmu = mu; 
    upcovar = covar;
else option == 2
    mu_min_j = n/(n-1)*(mu - 1/n*delrow);
    upmu = mu - 1/n*data(delrow,:) + 1/n*addrow;
    upcovar = covar - 1/n*((mu_min_j - data(delrow,:))'*(mu_min_j - data(delrow,:))) + 1/n*((mu_min_j - addrow)'*(mu_min_j - addrow))
end


