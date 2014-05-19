function result=qnm(x)

%QNM is a scale estimator which does not require an auxiliary location estimate.
% Essentially it is the first quartile of all pairwise distances between
% two data points. Its definition is given by 
%      qn(x)= c_n 2.2219{|x_i - x_j|; i<j}_(k)
% with c_n a small sample correction factor to make the qn unbiased at the
% normal distribution. It can resist 50% outliers.
% If x is a matrix, the scale estimate is computed on the columns of x. The
% result is then a row vector. If x is a row or a column vector,
% the output is a scalar.
% 
% Reference:
%   Rousseeuw P.J., Croux C. (1993), 
%   "Alternatives to the median absolute deviation", 
%   Journal of the American Statistical Association, 88, 1273-1283. 
%
% Required input argument:
%    x: either a data matrix with n observations in rows, p variables in columns
%       or a vector of length n.
% 
% I/O: result=qnm(x);
%
% This function is fully written in Matlab code. For a faster computation, use
% the qn.m function which calls qn.dll.
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by S. Verboven and Symen De Jong.
% Last Update: 03 April 2006

[n,p] = size(x);
if n==1 & p~=1
    x=x';
    [n,p]=size(x);
end
h = floor(n/2)+1;
Qn= zeros(1,p);
for j=1:p 
   dist = abs(repmat(x(:,j),1,n)-repmat(x(:,j)',n,1));
   d = sort(dist(triu(ones(n),1)>0));
   Qn(j) = 2.2219*d(h*(h-1)/2);
end
if n<=9
   switch n
   case 2
      dn=0.399;
   case 3
      dn=0.994;
   case 4
      dn=0.512;
   case 5
      dn=0.844;
   case 6
      dn=0.611;
   case 7
      dn=0.857;
   case 8
      dn=0.669;
   case 9
      dn=0.872;
   end
else
  if mod(n,2)==1
     dn=n/(n+1.4);
  end
  if mod(n,2)==0
     dn=n/(n+3.8);
  end
end
result=dn*Qn;     