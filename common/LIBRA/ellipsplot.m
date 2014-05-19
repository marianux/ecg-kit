function ellipsplot(center,covar,data,dist,nid,labx,laby)

%ELLIPSPLOT plots the 97.5% tolerance ellipse of the bivariate data set
% (data). The ellipse is defined by those data points whose distance (dist) 
% is equal to the squareroot of the 97.5% chisquare quantile with 2 degrees of
% freedom. 
%  
% Required input arguments:
%   center  : estimate of the center of the data set
%    covar  : estimate of the covariance matrix of the data set 
%     data  : the two-dimensional data matrix
%     dist  : distance needed to flag data points outside the ellipse
%
% Optional input arguments:
%      nid  : number of data points with largest distance
%             to be identified (default value: 3)
%     labx  : a label for the x axis (default value: 'X1')
%     laby  : a label for the y axis (default value: 'X2')
%
% I/O: ellipsplot(center,covar,data,dist,nid,labx,laby)
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Last update: 23/10/2003

set(gcf,'Name', '97.5% Tolerance ellipse', 'NumberTitle', 'off');
n=length(dist);
if size(data,2)~=2
     disp('The tolerance ellipse is only drawn for two-dimensional data sets')
else
    if nargin==4
        nid=3;
        labx='X1';
        laby='X2';
    elseif nargin==5
        labx='X1';
        laby='X2';
    elseif nargin==6
        laby='X2';
    end
     deter=covar(1,1)*covar(2,2)-covar(1,2)^2;
     ylimit=sqrt(7.37776*covar(2,2));
     y=-ylimit:0.005*ylimit:ylimit;
     sqtdi=sqrt(deter*(ylimit^2-y.^2))/covar(2,2);
     sqtdi([1,end])=0;
     b=center(1)+covar(1,2)/covar(2,2)*y;
     x1=b-sqtdi;
     x2=b+sqtdi;
     y=center(2)+y;
     ellip=[x1,x2([end:-1:1]);y,y([end:-1:1])]';
     xmin=min([data(:,1);ellip(:,1)]);
     xmax=max([data(:,1);ellip(:,1)]);
     ymin=min([data(:,2);ellip(:,2)]);
     ymax=max([data(:,2);ellip(:,2)]);
     xmarg=0.05*abs(xmax-xmin);
     ymarg=0.05*abs(ymax-ymin);
     xmin=xmin-xmarg;
     xmax=xmax+xmarg;
     ymin=ymin-ymarg;
     ymax=ymax+ymarg;
     x1=data(:,1)';
     x2=data(:,2)';
     plot(x1,x2,'o')
     xlabel(labx)
     ylabel(laby)
     xlim([xmin,xmax]);
     ylim([ymin,ymax]);
     box on
     [dist,ind]=sort(dist);
     ind=ind(n-nid+1:n)';
     text(x1(ind),x2(ind),int2str(ind));   
     title('Tolerance ellipse (97.5%)');         
     line(ellip(:,1),ellip(:,2),'Color','r');
 end