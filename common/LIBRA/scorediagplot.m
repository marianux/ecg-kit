 function scorediagplot(Xdist,Odist,k,cutoffx,cutoffo,attrib,labsd,labod)
 
%SCOREDIAGPLOT plots the score outlier map. The score distances (SD) and 
% orthogonal distances (OD) calculated in a PCA analysis are plotted 
% on the x and y-axis, respectively. 
% The cutoff values marked by a red line indicate which observations are
% outlying with respect to the majority of the data. The terminology used is 
% as indicated in the tabular:
%                   small SD          | large SD
%      ----------------------------------------------------------
%      small OD | regular point       | good PCA-leverage point
%      large OD | orthogonal outlier  | bad PCA-leverage point
%
% For more details see: 
%    Hubert, M., Verboven, S. (2003),
%    "A robust PCR method for high-dimensional regressors",
%    Journal of Chemometrics, 17, 438-452. 
%
% Required input arguments:
%  Xdist    : Score distances out of a PCA analysis
%  Odist    : Orthogonal distances out of a PCA analysis
%  k        : Number of principal components used in a PCA/PCR analysis
%  cutoffx  : Cutoff value for the score distances 
%             (=97.5% quantile of a chisquare with k degrees of freedom)
%  attrib   : string identifying the used method =  RPCR, ROBPCA, CPCR,
%             CPCA, MCDREG...
%
% Optional inputs
%   labsd   : number of displayed points with largest score distance (default = 3)
%   labod   : number of displayed points with largest orthogonal distance (default=3)
%
% I/O: scorediagplot(out.sd,out.od,out.k,out.cutoff.sd,out.cutoff.od,'RPCR',3,5)
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by S. Verboven
% Last update: 20/06/2003

%Checking input
if nargin==7
    labod=3;
end
if nargin==6
    labsd=3;
    labod=3;
end
if nargin<6
    error('no class identifier is given.')
end

%INITIALISATION%
%in case k=rank(T)=r then OD = zero 
NihilO=(Odist <= 1.e-06); 
Odist(NihilO)=0;
NihilX=(Xdist <= 1.e-06);
Xdist(NihilX)=0;
if ~any(Odist)     %OD=0: x=index and y=Xdist
    x=(1:length(Xdist))';
    y=Xdist;
    difflabel=1; %change name of the axis-labels
    if ~any(Xdist)
        axes
        box on
        text(0.05,0.5,'The score outliermap is not drawn when all OD and SD are zero.','color','r')
        return
    end
elseif ~any(Xdist) % SD=0: x=index and y=Odist
    x=(1:length(Xdist))';
    y=Odist;
    difflabel=2; %change name of the axis-labels
else
    x=Xdist;
    y=Odist; 
    difflabel=0;
end
quantx=cutoffx;
set(gcf,'Name', 'Score outlier map', 'NumberTitle', 'off');
hold on
xmin=0;
xmax=max([x; cutoffx]);
ymin=0;
ymax=max([y; cutoffo]);
xmarg=0.06*(xmax-xmin);
ymarg=0.06*(ymax-ymin);
xmin=xmin-xmarg;
xmax=xmax+xmarg;
ymin=ymin-ymarg;
ymax=ymax+ymarg;
n=length(x);
plot(x,y,'o');
if difflabel==1 %in case OD = 0
    ymax=max(ymax,cutoffx+ymarg);
    xlabel('Index')
    ylabel(['Score distance (',num2str(k),' LV)'])
    line('Xdata',[xmin xmax] ,'Ydata',[cutoffx cutoffx],'Linestyle','-','Color','r');
    xlim([xmin,xmax]);
    ylim([ymin,ymax]);
    box on
    plotnumbers(x',y,0,labsd,2);
    title(attrib)
    hold off
elseif difflabel==2 %case SD = 0
    xmax=max(xmax,cutoffo+xmarg);
    xlabel('Index')
    ylabel(['Orthogonal distance (',num2str(k),' LV)'])
    line('Xdata',[xmin xmax] ,'Ydata',[cutoffo cutoffo],'Linestyle','-','Color','r');
    xlim([xmin,xmax]);
    ylim([ymin,ymax]);
    box on
    plotnumbers(x',y,0,labod,2);
    title(attrib)
    hold off
else
    xlabel(['Score distance (',num2str(k),' LV)']);
    ylabel('Orthogonal distance');
    line([cutoffx cutoffx],[ymin ymax] ,'Linestyle','-','Color','r');
    line('Xdata',[xmin xmax] ,'Ydata',[cutoffo cutoffo], 'Linestyle','-','Color','r');
    xlim([xmin,xmax]);
    ylim([ymin,ymax]);
    box on
    plotnumbers(x',y,labsd,labod,2);
    title(attrib)
    hold off
end
