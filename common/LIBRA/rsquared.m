function [R2,out]=rsquared(x,y,kmax,method,h)

%RSQUARED calculates the R-squared value of the robust or classical PCR/PLS analysis. This function  
% is used in rpcr.m and rsimpls.m.
%
% The Robust R-squared is described in:
%    Hubert, M., Verboven, S. (2003),
%    "A robust PCR method for high-dimensional regressors",
%    Journal of Chemometrics, 17, 438-452.
%
% The required input arguments
%    x    : the explanatory variables
%    y    : the response variables
%    kmax : the number of components to choose in the ROBPCA method 
%  method : the method for which the R2 has to be computed. It can be 'RSIMPLS' or 'RPCR'. 
%
%Optional input arguments
%    h    : the quantile used in RPCR/RSIMPLS
%
% I/O: R2=rsquared(x,y,10)
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by S.Verboven on 31-10-2002
% Last updated on 18-02-2003


if nargin<4
    error('Missing input arguments')
end
[n,p]=size(x);
[n,q]=size(y);
cutoffWeights = sqrt(chi2inv(0.975,q));

if nargin<5
    h=floor(0.75*n);
end

count=1;
if q>1
    while count*q+q+(q*(q+1)/2) <= h
        count=count+1;
    end
else
    while count+2<h
        count=count+1;
    end
end
ktot=max(1,min(count-1,kmax));
if count-1<kmax
    disp(['Warning (rpcr): The maximal number of components is set to ', num2str(ktot),' to avoid overfitting.'])
end
weight=zeros(n,1);
bound=floor((n+ktot+q)/2);
h=max(h,bound);
switch method
case 'RPCR'
    outcv = cvRpcr(x,y,ktot,0,h);
case 'RSIMPLS'
    outcv = cvRsimpls(x,y,ktot,0,h);
end
R2 = outcv.R2;
W = outcv.outWeights.w_min;
rss = outcv.rss;

% place the two figures next to each other in the middle of the screen
bdwidth=5;
topbdwidth=30;
set(0,'Units','pixels');
scnsize=get(0,'ScreenSize');
pos1=[bdwidth, 1/3*scnsize(4)+bdwidth, scnsize(3)/2-2*bdwidth, scnsize(4)/2-(topbdwidth+bdwidth)];
pos2=[pos1(1)+scnsize(3)/2, pos1(2), pos1(3), pos1(4)];

figure('Position',pos1)
set(gcf,'Name', 'R-squared curve', 'NumberTitle', 'off');
plot(1:ktot,R2,'o-');
xlabel('Number of components');
ylabel('R2');
title(method)
ylabel('Robust R^2-value');
set(gca,'XTick',1:1:ktot)

figure('Position',pos2)
set(gcf,'Name', 'Square Root of Residual Sum of Squares curve', 'NumberTitle', 'off');
plot(1:ktot,sqrt(rss),'o-');
xlabel('Number of components');
title(method)
ylabel('Square Root of Robust RSS-value');
set(gca,'XTick',1:1:ktot)

k=input('How many components would you like to retain? ');
out.k=k;
out.weights=W(:,1:k);
out.rss = rss;

% closing the figures.
% close
% close

    
    
    
    
        