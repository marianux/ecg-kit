function residualplot(x,residuals,attrib,labx,laby,nid)

%RESIDUALPLOT plots the residuals from a regression analysis versus x
%
% Required input arguments:
%          x : the vector to be plotted on de x-axis
%  residuals : the residuals
%     attrib : string identifying the used method =  'LS', 'MLR', 'LTS', 'MCDREG' 
%
% Optional input arguments:
%       labx : a label for the x-axis (default: ' ')
%       laby : a label for the y-axis (default:' ')
%        nid : number of points to be identified in plots (default: 3)
%
% I/O: residualplot(x,residuals,attrib,labx,laby,nid)
%
% Last update: 07/04/2004
%

set(gcf,'Name', 'Residual plot', 'NumberTitle', 'off');
n=length(residuals);
if nargin<2
    error('A required input argument is missing.')
elseif nargin==3
    labx='';
    laby='';
    nid=3;
elseif nargin==4
    laby='';
    nid=3;
elseif nargin==5
    nid=3;
end

plot(x,residuals,'bo')
xlabel(labx);
ylabel(laby);
ord=abs(residuals);
[ord,ind]=sort(ord);
ind=ind(n-nid+1:n);
text(x(ind),residuals(ind),int2str(ind));
v=axis;
quant=sqrt(chi2inv(0.975,1));
line([v(1),v(2)],[quant,quant],'color','r');
line([v(1),v(2)],[0,0],'color','r');
line([v(1),v(2)],[-quant,-quant],'color','r');
xlim=([min(x), max(x)]);
ymin=min([-3 min(residuals)+0.05*min(residuals)]);
ymax= max([3 max(residuals)+0.05*max(residuals)]);
ylim([ymin,ymax]);
title(attrib)
box on