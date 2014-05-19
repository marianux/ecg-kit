function ddplot(x,y,cutoff,attrib,nid)

%DDPLOT is the distance-distance plot as introduced by Rousseeuw and Van
% Zomeren (1990, JASA, 85, 633-639). The Robust distances based on the MCD (mcdcov.m) 
% are plotted against the Mahalanobis distances. Cutoff lines permit the
% classification of outliers.
%
% Required input arguments:
%         x   : a vector containing the mahalanobis distances
%         y   : a vector containing the robust distances
%    cutoff   : the cutoff value for the distances
%
% Optional input arguments:
%       nid   : number of points to be identified in plots
%               (Default value: 3)
%
% I/O: ddplot(x,y,cutoff,nid)
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Last update: 24/11/2003

set(gcf,'Name', 'Distance-distance plot', 'NumberTitle', 'off');
if nargin==3
    nid=3;
end
ymax=max([max(y),cutoff,2.5])*1.05;
xmax=max([max(x),cutoff,2.5])*1.05;
plot(x,y,'o')
xlabel('Mahalanobis distance');
ylabel('Robust distance');
title(attrib)
xlim([-0.01*xmax,xmax]);
ylim([-0.01*ymax,ymax]);
box on
plotnumbers(x,y,0,nid,1);
line(repmat(max([cutoff,2.5]),1,2),[-0.01*ymax,ymax],'Color','r');
line([-0.01*xmax,xmax],repmat(max([cutoff,2.5]),1,2),'Color','r');
hold on
plot([-0.01*xmax,min([xmax,ymax])],[-0.01*ymax,min([xmax,ymax])],':','Color','g');
hold off
            
            
            
