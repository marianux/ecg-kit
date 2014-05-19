function normqqplot(y,class)

%NORMQQPLOT produces a Quantile-Quantile plot in which the vector y is plotted against 
% the quantiles of a standard normal distribution.
%
% Required input arguments:
%       y : row or column vector
%
% Optional input arguments:
%    class : a string used for the y-label and the title(default: ' ')
%
% This function is part of LIBRA: the Matlab library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
%Written by Nele Smets
% Last update: 12/03/2004

set(gcf,'Name', 'Normal QQ-plot', 'NumberTitle', 'off');
if nargin==1
    class=' ';
end
n=length(y);
for i=1:n
    normalquantile(i)=norminv((i-1/3)/(n+1/3),0,1);
end
y=sort(y);
plot(normalquantile,y,'bo')
xlabel('Quantiles of the standard normal distribution');
ylabel(['Standardized ',class,' residual']);
title(class);
box on
