function regresdiagplot(xdist,ydist,cutoffx,cutoffy,k,multi,attrib,labsd,labresd)

%REGRESDIAGPLOT makes a regression outlier map. 
% The score distances (SD)/mahalanobis distances (MD)/robust distances (RD) and 
% residual distances (ResD) calculated in a PCR or PLS analysis, or in a
% univariate or multivariate regression are plotted on the x and y-axis, respectively. 
% The cutoff values marked by a red line indicate which observations are
% outlying with respect to the majority of the data. The terminology used is 
% as indicated in the tabular:
%                         small SD/MD/RD  | large SD/MD/RD
%      ----------------------------------------------------------
%      small ResD/StdRes | regular point      | good leverage point
%      large ResD/StdRes | vertical outlier   | bad leverage point
%
% For more details see: 
%    Hubert, M., Verboven, S. (2003),
%    "A robust PCR method for high-dimensional regressors," 
%    Journal of Chemometrics, 17, 438-452. 
%
% Required input arguments:
%
%  xdist   : Score distances from a PCA analysis or
%            Mahalanobis distances from LS or MLR analysis or
%            Robust distances from LTS or MCDREG analysis, 
%  ydist   : Residual distances of the PCR method or
%            Standardized residuals of one of the regression methods,
%  cutoffx : Cutoff value for the SD/MD/RD,
%  cutoffy : Cutoff value for the residual-distance,
%  k       : Number of  principal components used in the PCR/PLS method, or zero
%            if not available.
%  multi   : 0 for univariate analysis, 1 for multivariate analysis
%  attrib  : String identifying the method used =  'LS', 'MLR', 'LTS', 'MCDREG', 'RPCR', 
%           'CPCR', 'RSIMPLS', 'CSIMPLS'
%
% Optional inputs:
%   labsd    : number of displayed points with largest distance on x-axes (default = 3)
%   labresd  : number of displayed points with largest distance on y-axes (default = 3)
%   
%
% I/O: regresdiagplot(out.sd,out.rd,out.cutoff.sd,out.cutoff.rd,k,multi,attrib,labsd,labresd)
%
% Example: regresdiagplot(out.sd,out.rd,out.cutoff.sd,out.cutoff.rd,out.k,1,'RPCR',5,5)
%          regresdiagplot(out.md,out.stdres,out.cutoffmd,0,0,0,'LS',5,4) 
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Sabine Verboven 
% Created 25 February 2002
% Revised : 12/02/2004
%

if nargin==8
    labresd=3;
end
if nargin==7 
    labsd=3;
    labresd=3;
end
if nargin<7
    error('A required input variable is missing!')
end
    
if ischar(xdist) | ischar(ydist)
    mess=sprintf(['Warning: A singularity was detected in the analysis. ',...
            '\n No regression outlier map available.']);
    disp(mess)
    return
end
%all LTS-analysis in RPCR are intercept included!!!
if multi==1 %multivariate analysis
    x=xdist; %distances on x-axis: SD, MD or RD
    y=ydist; %residual distances
    quanty=cutoffy;
    quantx=cutoffx;
else  %univariate analysis
    quantx=cutoffx; 
    quanty=sqrt(chi2inv(0.975,1)); 
    x=xdist;   %distances on x-axis: SD, MD, RD
    y=ydist;   %standardized residuals!!!
end

plot(x,y,'o');
set(gcf,'Name', 'Regression outlier map', 'NumberTitle', 'off');
hold on 
xmin=0;%max(0,min(x));%-0.5;
xmax=max([quantx max(x)])+0.5;

if multi==1
    ymin=0;%max(0,min(y));%-1;
    if size(y,2)==1
        y=y';
    end
    ymax=max([y,quanty])+1;
    xlim([xmin,xmax]);
    ylim([ymin,ymax]);
    ylimi=get(gca,'Ylim');
    line([quantx,quantx],[ylimi(1),ylimi(2)],'Linestyle','-','Color','r');
    xlimi=get(gca,'Xlim');
    line([xlimi(1),xlimi(2)],[quanty,quanty],'Linestyle','-','Color','r');
    if labsd
        plotnumbers(x',y',labsd,labresd,2)%sort by  x and  y
    end
    ylabel('Residual distance')
else
    ymin=min([-4 min(y)])-0.5;
    ymax=max([4 max(y)])+0.5;
    xlim([xmin,xmax]);
    ylim([ymin,ymax]);
    ylimi=get(gca,'Ylim');
    ylimi=[min(ylimi(1),min([-4 min(y)]))-0.5,max(ylimi(2),max([4 max(y)]))+0.5];
    line([quantx,quantx],[ylimi(1),ylimi(2)],'Linestyle','-','Color','r');
    xlimi=get(gca,'Xlim');
    xlimi=[min(xlimi(1),min(x))-1,max(xlimi(2),max([quantx max(x)]))+1];
    line([xlimi(1),xlimi(2)],[quanty,quanty],'Linestyle','-','Color','r');
    line([xlimi(1),xlimi(2)],[-quanty,-quanty],'Linestyle','-','Color','r');
    if labsd
        plotnumbers(x',y',labsd,labresd,3)%sort by x and abs(y)
    end
    ylabel('Standardized residual')
end

box on

switch attrib
    case 'LS'
        xlabel('Mahalanobis distance');
        ylabel('Standardized LS residual');
    case 'MLR'
        xlabel('Mahalanobis distance');
    case 'LTS'
        xlabel('Robust distance computed by MCD');
        ylabel('Standardized LTS residual');
    case 'MCDREG'
        xlabel('Robust distance computed by MCD');
    case {'CPCR'}
        xlabel(['Score distance (',num2str(k),' LV)']);
    case {'RPCR'}
        xlabel(['Score distance (',num2str(k),' LV)']);
    case {'CSIMPLS'}
        xlabel(['Score distance (',num2str(k),' LV)']);
    case {'RSIMPLS'}
        xlabel(['Score distance (',num2str(k),' LV)']);
     
end
title(attrib)       
hold off
