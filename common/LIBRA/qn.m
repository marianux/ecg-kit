function result=qn(x)

%QN is a scale estimator which does not require an auxiliary location estimate.
% Essentially it is the first quartile of all pairwise distances between
% two data points. Its definition is given by 
%      qn(x)= c_n 2.2219{|x_i - x_j|; i<j}_(k)
% with c_n a small sample correction factor to make the qn unbiased at the
% normal distribution. It can resist 50% outliers.
% If x is a matrix, the scale estimate is computed on the columns of x. The
% result is then a row vector. If x is a column vector, the output is a scalar.
% 
% Reference:
%   Rousseeuw P.J., Croux C. (1993), 
%   "Alternatives to the median absolute deviation", 
%   Journal of the American Statistical Association, 88, 1273-1283. 
%
% Required input argument:
%    x: either a data matrix with n observations in rows, p variables in columns
%       or a column vector of length n.
% 
% I/O: result=qn(x);
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% This is help on the C MEX-file "qn.dll". If full Matlab code is required, 
% use the qnm.m function.
%
% Original Fortran Code by C. Croux, translated to MATLAB by S. Verboven.
% Last update : 13/07/2000
