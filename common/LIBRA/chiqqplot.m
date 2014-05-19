function chiqqplot(y,p,class)

%CHIQQPLOT produces a Quantile-Quantile-plot of the vector y 
% versus the square root of the quantiles of the chi-squared distribution.
%
% Required input arguments:
%       y  : row or column vector 
%       p  : degrees of freedom of the chi-squared distribution 
%
% Optional input argument:
%    class : a string used for the y-label and the title(default: ' ')
%
% I/O: chiqqplot(y,p,class)
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Last update: 23/10/2003

set(gcf,'Name', 'Chisquare QQ-plot', 'NumberTitle', 'off')
n=length(y);
if nargin==2
    class='';
end
for i=1:n
     x(i)=chi2inv((i-1/3)/(n+1/3),p);
end
x=sqrt(x);
y=sort(y);
plot(x,y,'o')
xlabel('Square root of the quantiles of the chi-squared distribution');
if strcmp(class,'MCDCOV')
    ylabel('Robust distance');
elseif strcmp(class,'COV')
    ylabel('Mahalanobis distance');
else
    ylabel('Distance');
end
title(class);
