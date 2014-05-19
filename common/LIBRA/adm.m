function result=adm(x)

%ADM is a scale estimator given by the Average Distance to the Median.
% It is defined as
%    adm(x) = ave(|x_i - med(x)|)
% If x is a matrix, the scale estimate is computed on the columns of x. The
% result is then a row vector. If x is a row or a column vector,
% the output is a scalar.
%
% The ADM is also described in:
%   Rousseeuw, P.J. and Verboven, S. (2002),
%   "Robust Estimation in Very Small Samples", 
%   Computational Statistics and Data Analysis, 40, 741-758.
%
% Required input argument:
%    x : either a data matrix with n observations in rows, p variables in columns
%        or a vector of length n.
%
% I/O: out=adm(x);
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Sabine Verboven
% Last revision 28/08/03 by Nele Smets

[n,p]=size(x); 

if n==1 & p==1
    result=0;  %if X is of size 1x1, the scale estimate is equal to 0
    return
elseif  n==1      
    x=x';   %if X is row vector, transpose to a column vector
    n=p;
    p=1;
end

for i=1:p
    X=x(:,i);
    result(i)=mean(abs(X-median(X)));
end

