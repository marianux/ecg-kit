function result=rrmse(x,y,h,kmax,attrib,plots,k,weight,res)

%RRMSE calculates the robust RMSECV and/or the robust RMSEP-value 
% for RPCR and RSIMPLS.
%
% The robust RMSECV is described in:
%  
%  Engelen, S., Hubert, M. (2005),
%  "Fast model selection for robust calibration methods",
%  Analytica Chimica Acta, 544, 219-228.
%
% Required input arguments:
%        x : the regressors
%        y : the response variables
%        h : the quantile used in RPCR and RSIMPLS.
%     kmax : the maximal number of components to be used.
%   attrib : the name of the analysis 'RSIMPLS', or 'RPCR'  
%
% Optional input arguments:
%    plots : 0/1 = no plot / plot (default) of the RMSECV values 
%        k : the optimal number of components chosen by cross validation.
%            if k is different from zero, the RMSEP will be calculated
%            (default k=0)
%   weight : only needed when calculating RMSEP value
%    res   : the residuals of each left-out observation for every k.
%
% I/O: result=rrmse(x,y,h,kmax,attrib,plots,k,weight,res);
%
% Example: h=0.75*size(x,1);
%          result=rrmse(x,y,h,10,'RPCR'); 
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Created by Karlien Vanden Branden on 05-07-2002
% Revisions by Sabine Verboven, Sanne Engelen
% Last update on 28-04-2004

%%%%%%%%% Initialisation %%%%%%%%%%%%%%
if nargin>=7
    in.pred=1;
else 
    in=struct('');
end
if nargin<7 
    k=0; 
end
if nargin<6
    k=0;
    plots=1;
end
if nargin<5
    error('Missing one or more input variables.')
end
if k<0 %preventing negative input of chosen number of PC's
    k=0;
end
   
[n,p]=size(x);
[n,q]=size(y);

%%%%%Defining the maximum number of principal components to calculate
count=1;
if q>1
    while (count*q+q+(q*(q+1)/2))<=h
        count=count+1;
    end
else
    while (count+2)<=h
        count=count+1;
    end
end
ktot=max(1,min(count-1,kmax));
%%%%%%%%%%%%%%%% MAIN PART %%%%%%%%%%%%%%%%%%%%%%%%
if (k~=0) & nargin>7
    disp(['The RMSEP value is based on ', num2str(k),' scores.'])
    weight=weight(:,k);
    res=res(:,(k-1)*q + 1:k*q).^2;
    res=repmat(weight,1,q).*res;
    vRMSECV2=sum(res,1)/sum(weight); %1xn -> 1x1 
    out.rmsep=sqrt(sum(vRMSECV2)/q);  
else
    % Cross-validation 
    if k==0
        disp(['Cross-validation is now performed.'])
    else
        disp(['The RMSEP value is based on ', num2str(k),' components.'])
    end
    
    resCV = crossvalid(attrib,x,y,ktot,h,k);
    out.R2 = resCV.R2;
    out.rss = resCV.rss;
    
    if plots & (k==0)
        disp(['The robust RMSECV-values: ', num2str(resCV.rmsecv)])
        figure
        set(gcf,'Name', 'Robust Component Selection plot', 'NumberTitle', 'off');
        plot(1:ktot,resCV.rmsecv,'o-');
        hold on
        plot(1:ktot,sqrt(0.5*resCV.rmsecv.^2 + 0.5*resCV.rss),'r*--')
        plot(1:ktot,sqrt(resCV.rss),'g>-')
        xlabel('Number of components');
        set(gca,'XTick',1:1:ktot)
        ylabel('RCS value');
        title(attrib)
        legend('\gamma = 1 (CV)','\gamma = 0.5','\gamma = 0 (RSS)')
        hold off
        kout=input(['How many components would you like to retain? ']);
        out.k=kout;
    else
        out.k=k;
    end
    if k==0 %output needed for the calculation of rmsep
        out.weight=resCV.outWeights.weightsk;
        out.res=resCV.residu;
        out.rmsecv=resCV.rmsecv;
    else
        out.rmsep=resCV.rmsep;
        out.k=k;
    end
end

result = out;
%%%%%%%%%%%%%%%%%%%%%%%%%subfunction%%%%%%%%%%%%%%%%%%%%%%

function out=crossvalid(attrib,x,y,kmax,h,k)

switch attrib
case 'RPCR'
    out = cvRpcr(x,y,kmax,1,h,k);   
case 'RSIMPLS'
    out = cvRsimpls(x,y,kmax,1,h,k);   
end     


