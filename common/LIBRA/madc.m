function result=madc(x)

%MADC is a scale estimator given by the Median Absolute Deviation 
% with finite sample correction factor.
% It is defined as 
%           mad(x)= b_n 1.4826 med(|x_i - med(x)|)
% with b_n a small sample correction factor to make the mad unbiased at the
% normal distribution. It can resist 50% outliers. 
% If x is a matrix, the scale estimate is computed on the columns of x. The
% result is then a row vector. If x is a row or a column vector,
% the output is a scalar.
% 
% Required input argument:
%    x: either a data matrix with n observations in rows, p variables in columns
%       or a columnn vector of length n.
% 
% I/O: result=mad(x);
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by S.Verboven 
% Last revised: 22/12/03 

[n,p]=size(x);

if n==1 & p==1
    result=0;  %when X is a one by one matrix, all scale estimators must be equal to 0
    return
elseif n==1        
    x=x';   %we only want to work with column vectors
    n=p;
    p=1;
end
bn=fcorfac(x); %calculating finite sample correction
t_0=median(x);
result=median(abs(x - repmat(median(x),n,1)))*1.4826*bn; 
if ~all(result)
    result(result==0)=t_0(result==0);
end

%---------
function bn=fcorfac(Z)

[n,p]=size(Z);
switch n
case 2
   bn=1.196;
case 3
   bn=1.495;
case 4
   bn=1.363;
case 5
   bn=1.206;
case 6
   bn=1.200;
case 7
   bn=1.140;
case 8
   bn=1.129;
case 9
   bn=1.107;
end
if n>9
   bn=n/(n-0.8);
end












