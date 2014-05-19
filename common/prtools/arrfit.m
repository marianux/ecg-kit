function [beta,msr] = arrfit(X,y,lambda,precision)

% ---------------------------------------------------------------------------------
% [beta,msr] = arrfit(X,y,lambda,precision)
%
% Adaptive Ridge Regression linear fit to data 
% Finds the coefficients BETA, of the linear fit to the data, 
%    X(i,:)*BETA ~= y(i), 
% minimizing the following expression:
%   sum((X*BETA-y).^2) + lambda * sum(abs(BETA))^2
%
% INPUT:
%           X: (Nsamples X Nfeatures) vector or matrix of input data
%           y: (Nsample X 1) target values or output data
%      lambda: (scalar or vector) (default=1) penalty coefficients 
%              If lambda is a vector, each column of beta and msr
%              corresponds to the respective value of lambda. 
%   precision: (scalar) (default=1e-2) (optional) measure of the 
%              absolute and relative precisions required for beta.  
%
% OUTPUT:
%        beta: (Nfeatures x 1) the regression coefficients 
%         msr: the mean squares residuals.
%
% 22/06/98 Y. Grandvalet 
% @inproceedings{Grandvalet98a,
%    AUTHOR = "Grandvalet, Y.",
%     TITLE = "Least absolute shrinkage is equivalent to quadratic penalization ",
% BOOKTITLE = "ICANN 98",
%    EDITOR = "Niklasson, L. and Bod{\'e}n, M. and Ziemske, T.",
%    VOLUME = "1",
%     PAGES = "201--206",
% PUBLISHER = "Springer",
%    SERIES = "Perspectives in Neural Computing",
%      YEAR = "1998"}
% ---------------------------------------------------------------------------------

if nargin < 4;
   precision = 1e-2;
   if nargin < 3;
      lambda = 1;
      if nargin < 2;
        error('ARRFIT requires at least two input arguments.');
      end;
   end;
end;
precision = precision.^2;
 
% Check that matrix (X) and vector (y) have compatible dimensions
 
[n,d]   = size(X);
[ny,dy] = size(y);
if ny~=n, 
    error('The number of rows in y must equal the number of rows in X.'); 
end 
if dy ~= 1, 
    error('y must be a vector, not a matrix'); 
end

% Check that (lambda) has correct dimensions
 
[nl,dl] = size(lambda);
if dl ~= 1 & nl ~= 1, 
    error('lambda must be a scalar or vector.');
end 
[nl] = max([nl,dl]);

% Check that (precision) has correct dimensions
 
if length(precision) ~= 1, 
    error('precision must be a scalar.');
end 

% Initializations

beta = zeros(d,nl); 

XX = (X'*X);
Xy = (X'*y);

for i=1:nl
    if  lambda(i)==Inf;
    beta(:,i) = zeros(d,1);
    else;
        Lambda = lambda(i)*ones(d,1); 
    U = chol(XX + diag(Lambda));
        betanew = U\(U'\Xy);
        stop    = 0;
        while (~stop);
             betaold      = betanew;
         normbetaold  = abs(betaold)./mean(abs(betaold));
             ind          = find( normbetaold > precision );
             Lambda(ind)  = (d*lambda(i))./normbetaold(ind);   
         betanew      = zeros(d,1); 
         U            = chol(XX(ind,ind) + diag(Lambda(ind)));
             betanew(ind) = U\(U'\Xy(ind));
             stop         = max( abs(betaold-betanew)./(1+abs(betanew)) ) < precision;
        end
        beta(:,i) = betanew;
    end;
end;

if nargout > 1
   msr = sum( (X*beta - y(:,ones(nl,1))).^2 )/n;
end;
