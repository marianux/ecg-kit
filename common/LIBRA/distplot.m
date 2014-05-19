function distplot(y,cutoff,class,nid)

%DISTPLOT plots the vector y versus the index.
% At the height of the cutoff value a red vertical line is plotted.
%
% Required input arguments:
%       y  : the vector to be plotted
%  cutoff  : a cutoff value 
%
% Optional input arguments:
%   class  : the class of the y-vector 
%     nid  : number of points to be identified (Default value: 3)
%
% I/O: distplot(y,cutoff,class,nid)
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at:  
%              http://wis.kuleuven.be/stat/robust.html
%
% Last update: 24/11/2003

set(gcf,'Name', 'Index plot of the distances', 'NumberTitle', 'off');
n=length(y);
if nargin==2
    laby='Distance';
    nid=3;
elseif nargin==3
    nid=3;
end
ymax=max([max(y),cutoff,2.5])*1.05;
plot(1:n,y,'o')
box on
xlabel('Index')
if strcmp(class,'MCDCOV')
    laby='Robust distance';
elseif strcmp(class,'COV')
    laby='Mahalanobis distance';
end
ylabel(laby)
xlim([-0.025*n,n*1.05]);
ylim([-0.025*ymax,ymax]);
plotnumbers(1:n,y,0,nid,1)
line([-0.025*n,n*1.05],repmat(max([cutoff,2.5]),1,2),'Color','r');
title([class]);