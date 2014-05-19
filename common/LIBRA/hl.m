function result=hl(x)

%HL computes the Hodges-Lehmann location estimate on the columns of x.
% The Hodges-Lehmann estimator is defined as 
%           hl(x)=med {(x_i + x_j)/2} with 1<=i<j<=n 
% It can resist about 29% outliers.
% If x is a matrix, the location estimate is computed on the columns of x. The
% result is then a row vector. If x is a row or a column vector,
% the output is a scalar.
%
% The behavior of the HL estimator in small samples is discussed in:
%   Rousseeuw, P.J. and Verboven, S. (2002),
%   "Robust estimation in very small samples", 
%   Computational Statistics and Data Analysis, 40, 741-758.
%
% Required input argument:
%    x: either a data matrix with n observations in rows, p variables in columns
%       or a vector of length n.
%
% I/O:  result=hl(x);
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by S. Verboven
% Last revision: 28/08/03 by N. Smets

[n,p]=size(x);

if n==1 & p==1
    result=x;     %when X is a one by one matrix, all location estimators must be equal to that matrix
    return
elseif n==1
    x=x';      %we only want to work with column vectors
    n=p;
    p=1;
end

if n==2
    result=mean(x,1);  % all location estimators must equal to the average for n=2
    return
end

for k=1:p
    m=0;
    y=zeros(1,n*(n-1)/2); %initializing help vector 
    X=x(:,k);
    for i=1:n              %calculating all possible pairwise means
        for j=(i+1):n
            m=m+1;
            y(m)=(X(i)+X(j))/2;
        end
    end
    result(k)=median(y);
    y=[];                 %clear help vector
end
  
