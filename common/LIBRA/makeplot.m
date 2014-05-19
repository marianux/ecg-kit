function makeplot(out,varargin)

%MAKEPLOT makes plots for the main functions. These figures can also be obtained 
% by setting 'plots = 1' in those functions.
%
% Required input:
%  out = a structure containing the output of one of the following classes: 
%        MCDCOV, LS, LTS, MLR, MCDREG, CPCA,CPCR, CSIMPLS, RAPCA, ROBPCA, 
%        RPCR, RSIMPLS, CDA, RDA
%
% Optional input: 
%  nameplot  :  0         : menu of plots     (default = 0) 
%              'all'      : all possible plots
%              'scree'    : scree plot
%              'pcadiag'  : Score outlier map
%              'regdiag'  : Regression outlier map
%              '3ddiag'   : 3D Regression outlier map
%              'robdist'  : A plot of the robust distances
%              'qqmcd'    : A qq-plot of the robust distances versus the quantiles of the chi-squared distribution
%              'dd'       : A DD-plot: Robust distances versus Mahalanobis distances
%              'ellipse'  : A tolerance ellipse
%              'resfit'   : Standardized LTS Residuals versus fitted values
%              'resindex' : Standardized LTS Residuals versus index
%              'qqlts'    : Normal QQ-plot of the LTS residuals 
%              'scatter'  : Scatterplot with LTS line 
%              'da'       : Tolerance ellipses of a discrimant analysis
%              'simca'    : Scatterplot with boundaries defined by the number of PC's from a simca method
%  labod   : number of points to be identified in score plots 
%            with largest orthogonal distance.
%  labsd   : number of points to be identified in score and regression plots 
%            with largest score distance.
%  labresd : number of points to be identified in regression plots 
%            with largest residual distance.
%  labmcd  : number of points to be identified in MCD plots
%            with the largest robust distance.
%  lablts  : number of points to be identified in LTS plots 
%            with the largest absolute standardized residual.
%  classic : In case of a robust analysis, classic can be set to 0 to avoid classical plots (default). 
%            If set to one, the inputargument 'out' must contain a field called 'classic' which is a structure. 
%
% I/O: makeplot(out,'nameplot',0,'labsd',3,'labod',3,'labresd',3,'classic',0)
%  The user should only give the input arguments that have to change their default value.
%  The name of the input arguments needs to be followed by their value. The order of 
%  the input arguments is of no importance.
%
% Example: outrpcr=rpcr(x,y,'plots',0,'classic',1); 
%          makeplot(outrpcr,'labsd',5, 'classic',1)
%          makeplot(outrpcr,'nameplot','scree')
%
%          outmcd=mcdcov(x);
%          makeplot(outmcd,'labmcd',6)
%          makeplot(outmcd,'labmcd',4,'nameplot','dd')
%
%          outls=ols(x,y);
%          makeplot(outls)
%
% This function is part of the Matlab Library for Robust Analysis (LIBRA),
% available at: 
%              www.wis.kuleuven.ac.be/stat/robust.html
%
% Written by Sabine Verboven  09/04/2003
% Last Updated : 29/01/2008
%
% Uses functions:  screeplot, scorediagplot, regresdiagplot,
%                  regresdiagplot3d, plotnumbers, plotnumbers3d,
%                  distplot, chiqqplot, ddplot, ellipsplot, residualplot, 
%                  normqqplot, lsscatter, menureg, menuscoreg,
%                  menupls, menucov, menuls,menuda, menusimca
%

%INTITIALISATION%
if isstruct(out)
    innames=fieldnames(out);
    if strmatch('class',innames,'exact')
        attrib=out.class;
    else
        error('The method had no class identifier.')    
    end
    if strmatch('classic',innames,'exact')
        classic=out.classic;
    elseif ismember(attrib,{'CPCA','CPCR','CSIMPLS','MLR','CDA','LS','CSIMCA'})
        classic= 1; 
    else 
        classic=0;
    end
else
    error('The first inputargument is not a structure.')
end
counter=1;
default=struct('nameplot',0,'labsd',3,'labod',3,'labresd',3,'labmcd',3,'lablts',3,'classic',1);
list=fieldnames(default);
options=default;
IN=length(list);
i=1;
%
if nargin>1
    %
    %placing inputfields in array of strings
    %
    for j=1:nargin-1
        if rem(j,2)~=0
            chklist{i}=varargin{j};
            i=i+1;
        end
    end 
    %
    %Checking which default parameters have to be changed
    % and keep them in the structure 'options'.
    %
    while counter<=IN 
        index=strmatch(list(counter,:),chklist,'exact');
        if ~isempty(index) %in case of similarity
            for j=1:nargin-2 %searching the index of the accompanying field
                if rem(j,2)~=0 %fieldnames are placed on odd index
                    if strcmp(chklist{index},varargin{j})
                        I=j;
                    end
                end
            end
            options=setfield(options,chklist{index},varargin{I+1});
            index=[];
        end
        counter=counter+1;
    end
    choice=options.nameplot;
    if choice~=0
        ask=0; 
    else       
        ask=1; %menu of plots 
    end
    labsd=options.labsd;
    labod=options.labod;
    labresd=options.labresd;
    labmcd=options.labmcd;
    lablts=options.lablts;
    classicplots=options.classic;
    if ~isstruct(classic) & options.classic==1
        mess=sprintf(['The classical output is not available. Only robust plots will be shown.\n',...
                'Please rerun the preceeding analysis with the option ''classic'' set to 1',...
              '\n if the classical plots are required.']);
        disp(mess)     
        classicplots=0;
    end
else
    ask=1; %menu of plots
    choice=0;
    labsd=3;
    labod=3;
    labresd=3;
    labmcd=3;
    lablts=3;    
    if ~isstruct(classic) & options.classic==1
       % mess=sprintf(['The classical output is not available. Only robust plots will be shown.\n',...
         %       'Please rerun the preceeding analysis with the option ''classic'' set to 1',...
           %     '\n if the classical plots are required.']);
        %disp(mess)     
        classicplots=0;
    else
        classicplots=1;
    end
end    

%to initialize the correct menu of plots
exitno=0;
switch attrib
case 'CPCA'
    exitno=4;
    if ismember({char(choice)},{'regdiag','3ddiag','robdist','qqmcd','dd','ellipse','resfit','resindex','qqlts','diag','scatter',...
                'da','simca'})
         error('That kind of plot is not available for this method.')
    end
case 'CPCR'
    exitno=6;
    if ismember({char(choice)},{'robdist','qqmcd','dd','ellipse','resfit','resindex','qqlts','diag','scatter','da','simca'})
         error('That kind of plot is not available for this method.')
    end   
case 'LS'
    exitno=7;    
     if ismember({char(choice)},{'scree','pcadiag','3ddiag','robdist','qqmcd','dd','ellipse','simca','da'})
         error('That kind of plot is not available for this method.')
    end
case 'MCDREG'
    exitno=3;
     if ismember({char(choice)},{'scree','pcadiag','3ddiag','robdist','qqmcd','dd','ellipse',...
             'resfit','resindex','qqlts','diag','scatter','da','simca'})
         error('That kind of plot is not available for this method.')
    end
case 'MLR'
    exitno=3;
     if ismember({char(choice)},{'scree','pcadiag','3ddiag','robdist','qqmcd','dd','ellipse',...
             'resfit','resindex','qqlts','diag','scatter','da','simca'})
         error('That kind of plot is not available for this method.')
    end
case 'LTS' 
    exitno=7;
    if ismember({char(choice)},{'scree','pcadiag','3ddiag',...
            'robdist','qqmcd','dd','ellipse','da','simca'})
         error('That kind of plot is not available for this method.')
    end
    if ~isfield(out,'rd')
        error('Please rerun the LTS-regression again with the option plots equal to 1.')
    end
case {'RAPCA','ROBPCA'}
    exitno=4;
   if ismember({char(choice)},{'regdiag','3ddiag','robdist','qqmcd','dd','ellipse','resfit',...
           'resindex','qqlts','diag','scatter','da','simca'})
        error('That kind of plot is not available for this method.')
   end
case 'RPCR' 
      exitno=6;
      if ismember({char(choice)},{'robdist','qqmcd','dd','ellipse','resfit','resindex','qqlts','diag','scatter',...
              'da','simca'})
         error('That kind of plot is not available for this method.')
    end   
case {'RSIMPLS','CSIMPLS'}
    exitno=5;
     if ismember({char(choice)},{'scree','robdist','qqmcd','dd','ellipse','resfit','resindex','qqlts','diag','scatter',...
             'da','simca'})
         error('That kind of plot is not available for this method.')
    end
case 'MCDCOV'
    if ismember({char(choice)},{'scree','pcadiag','regdiag','3ddiag','resfit','resindex','qqlts','diag','scatter',...
                'da','simca'})
         error('That kind of plot is not available for this method.')
    end
    exitno=6;
case {'RDA','CDA'}
    exitno=3;
    if ismember({char(choice)},{'scree','pcadiag','3ddiag','robdist','qqmcd','dd','ellipse',...
             'resfit','resindex','qqlts','diag','scatter','regdiag','simca'})
        error('That kind of plot is not available for this method.')
    end 
    if size(out.center,2)>2
        disp('Warning: Tolerance ellipses are only drawn for two-dimensional data sets.')
        return
    end
case {'CSIMCA','RSIMCA'}
    exitno=3;
    if ismember({char(choice)},{'scree','pcadiag','3ddiag','robdist','qqmcd','dd','ellipse',...
             'resfit','resindex','qqlts','diag','scatter','regdiag','da'})
        error('That kind of plot is not available for this method.')
    end 
    if size(out.pca{1}.P,1) > 3 
        disp('Warning: The dimension of the dataset is larger than 3.')
        return
    end
end
if exitno==0
    error(['Your attribute identifier must be one of the following names:',...
            'CPCA, RAPCA, ROBPCA, CPCR, RPCR, LS, LTS, MCDREG, RSIMPLS, CSIMPLS, MCDCOV, CDA,RDA,CSIMCA,RSIMCA'])
end

%plotting what is asked for 
if ask==0
    whichplot(out,choice,attrib,exitno,labsd,labod,labresd,labmcd,lablts,classic,classicplots)
else
    %make menu of plots
    while (choice~=exitno)
        switch attrib               
            case {'CPCA','ROBPCA','RAPCA'}
                choice=menuscore(out,attrib,exitno,labsd,labod,classic,classicplots);
            case {'MCDREG','MLR'}
                choice=menureg(out,attrib,exitno,labsd,labresd,classic,classicplots);
            case {'COV','MCDCOV'}
                choice=menucov(out,attrib,exitno,labmcd,classic,classicplots);
            case {'RSIMPLS','CSIMPLS'}
                choice=menupls(out,attrib,exitno,labsd,labod,labresd,classic,classicplots);
            case {'CPCR','RPCR'}
                choice=menuscoreg(out,attrib,exitno,labsd,labod,labresd,classic,classicplots);
            case {'LTS','LS'}
                 choice=menuls(out,attrib,exitno,labsd,labresd,lablts,classic,classicplots);    
            case{'CDA','RDA'}
                    choice=menuda(out,attrib,exitno,classic,classicplots);
            case{'CSIMCA','RSIMCA'}
                    choice=menusimca(out,attrib,exitno,classic,classicplots);
           end 
    end
end

%%%%%%%%%%%%% MAIN FUNCTION %%%%%%%%%%%%%%%
function whichplot(out,plotn,attrib,exitno,labsd,labod,labresd,labmcd,lablts,classic,classicplots)

%Initializing variables

switch attrib
    
case {'ROBPCA','RAPCA'}
    Xdist=out.sd;
    cutoffX=out.cutoff.sd;
    cutoffO=out.cutoff.od;
    OD=out.od;
    k=out.k;
    L=out.L;
    %classical output
    if isstruct(classic)
        Lcl=classic.L;
        Xdistcl=classic.sd;
        ODcl=classic.od;
        cutoffXcl=classic.cutoff.sd;
        cutoffOcl=classic.cutoff.od;
        attribcl='CPCA';
    end
case 'MCDREG'
    fitted=out.fitted;  
    Rdist=out.resd;
    cutoffR=out.cutoff.resd;
    if size(fitted,2)~=1 
        multi=1;
    else
        multi=0;
    end
    res=out.res;
    Xdist=out.rd;
    Rsquared=out.rsquared;
    cutoffX=out.cutoff.rd;
    Se=out.cov;
    k=size(out.slope,1);
    %classical output
    if isstruct(classic)
        attribcl=classic.class;
        fittedcl=classic.fitted;
        Rdistcl=classic.resd;
        standrescl=classic.stdres;
        cutoffRcl=classic.cutoff.resd;
        rescl=classic.res;
        Xdistcl=classic.md;
        cutoffXcl=classic.cutoff.md;
        Rsquaredcl=classic.rsquared;
        Secl=classic.cov;
    end
case 'MLR'
    fittedcl=out.fitted;
    Rdistcl=out.resd;
    cutoffRcl=out.cutoff.resd;
    if size(fittedcl,2)~=1
        multi=1;
    else
        multi=0;
    end    
    rescl=out.res;
    Xdistcl=out.md;
    cutoffXcl=out.cutoff.md;
    Rsquaredcl=out.rsquared;
    Secl=out.cov;
    k=size(out.slope,1);
    attribcl=out.class;
case 'RPCR'
    fitted=out.fitted; 
    Rdist=out.resd;
    cutoffR=out.cutoff.resd;
    if size(fitted,2)~=1
        multi=1;
    else
        multi=0;
    end
    res=out.res;
    Xdist=out.sd;
    Rsquared=out.rsquared;
    cutoffX=out.cutoff.sd;
    cutoffO=out.cutoff.od;
    OD=out.od;
    Se=out.cov;
    k=out.k;
    L=out.robpca.L;
    %classical output
    if isstruct(classic)
        attribcl='CPCR';
        Lcl=classic.cpca.L;
        Xdistcl=classic.sd;
        ODcl=classic.od;
        cutoffOcl=classic.cutoff.od;
        cutoffXcl=classic.cutoff.sd;
        Rdistcl=classic.resd; 
        cutoffRcl=classic.cutoff.resd;
    end
case 'LTS'
    resid=out.res;
    fitted=out.fitted;
    scale=out.scale;
    standres=resid/scale;
    n=length(resid);
    if isfield(out,'X')
        x=out.X;
        y=out.y;
    else
        x=0;
        y=0;
    end
    Xdist=out.rd;
    cutoffX=out.cutoff.rd;
    cutoffY=out.cutoff.rd;
    k=size(out.slope,1);
    %classical output
    if isstruct(classic)
        attribcl='LS';
        residcl=classic.res;
        fittedcl=classic.fitted;
        scalecl=classic.scale;
        standrescl=residcl/scalecl;
        Xdistcl=classic.md;
        cutoffXcl=classic.cutoff.md;
        ncl=length(residcl);
        if isfield(classic,'X')
            xcl=classic.X;
            ycl=classic.y;
        else
            xcl=0;
            ycl=0;
        end
    end
case 'LS'
    attribcl=attrib;
    residcl=out.res;
    fittedcl=out.fitted;
    scalecl=out.scale;
    Xdistcl=out.md;
    cutoffXcl=out.cutoff.md;
    standrescl=out.resd;
    ncl=length(residcl);
    if isfield(out,'X')
        xcl=out.X;
        ycl=out.y; 
    else
        xcl=0;
        ycl=0; 
    end
    k=size(out.slope,1);
case 'CPCA'
    Lcl=out.L;
    Xdistcl=out.sd;
    ODcl=out.od;
    cutoffOcl=out.cutoff.od;
    cutoffXcl=out.cutoff.sd;
    k=out.k;
    attribcl=attrib;
case 'CPCR'
    fitted=out.fitted;
    Lcl=out.cpca.L;
    k=out.k;
    Xdistcl=out.sd;
    ODcl=out.od;
    cutoffOcl=out.cutoff.od;
    cutoffXcl=out.cutoff.sd;
    Rdistcl=out.resd;
    cutoffRcl=out.cutoff.resd;
    if size(fitted,2)~=1
        multi=1; %multivariate analysis
    else
        multi=0; %univariate analysis
        stdrescl=out.resd;
        cutoffRcl=out.cutoff.resd;
    end   
    attribcl=attrib;
case 'RSIMPLS'
    fitted=out.fitted; 
    Rdist=out.resd;
    cutoffR=out.cutoff.resd;
    if size(fitted,2)~=1
        multi=1;
    else
        multi=0;
    end
    res=out.res;
    Xdist=out.sd;
    cutoffX=out.cutoff.sd;
    cutoffO=out.cutoff.od;
    OD=out.od;
    Se=out.cov;
    k=out.k;
    %classical output
    if isstruct(classic)
        attribcl='CSIMPLS';
        Xdistcl=classic.sd;
        ODcl=classic.od;
        cutoffOcl=classic.cutoff.od;
        cutoffXcl=classic.cutoff.sd;
        Rdistcl=classic.resd;
        cutoffRcl=classic.cutoff.resd;
    end
case 'CSIMPLS'
    fittedcl=out.fitted;
    if size(fittedcl,2)~=1
        multi=1;
    else
        multi=0;
    end
    Xdistcl=out.sd;
    ODcl=out.od;
    k=out.k;
    cutoffXcl=out.cutoff.sd;
    cutoffOcl=out.cutoff.od;
    attribcl=attrib;
    Rdistcl=out.resd;
    cutoffRcl=out.cutoff.resd;
case 'MCDCOV'
    if ~isempty(out.plane)
        disp('Warning (makeplot): The MCD covariance matrix is singular. No plots can be drawn.')
        return
    end
    covar=out.cov;
    md=out.md;
    rd=out.rd;
    if isfield(out,'X')
        data=out.X;
    else
        data=0;
    end
    center=out.center;
    cutoffR=out.cutoff.rd;  
    cutoffM=out.cutoff.md; 
    if isstruct(classic)
        if isfield(out,'X')
            datacl=out.X;
        else
            datacl=0;
        end
        centercl=out.classic.center;
        covcl=out.classic.cov; 
        mdcl=out.classic.md;
        attribcl=out.classic.class;
    end
case 'COV'
   %only available inside the mcdcov function.   
case 'CDA'
    if isfield(out,'x')
        xcl=out.x;
        groupcl=out.group;
    else
        xcl=0;
        groupcl=0;
    end
    methodcl=out.method;
    centercl=out.center;
    covcl=out.cov;
    classcl=out.class;
case 'RDA'
    if isfield(out,'x')
        x=out.x; group=out.group;
    else
        x=0;group=0;
    end
    method=out.method;
    center=out.center;
    covar=out.cov;
    class=out.class;
    if isstruct(classic)
        if isfield(out.classic,'x')
            xcl=out.classic.x;
            groupcl=out.classic.group;
        else
            xcl=0;
            groupcl=0;
        end
        methodcl=out.classic.method;
        centercl=out.classic.center;
        covcl=out.classic.cov;
        classcl=out.classic.class;
    end
case 'RSIMCA'
    if isstruct(classic)
        resultcl=out.classic;
    else
        resultcl=0;
    end
end

%determination of the number of plots (robust or/and classic) 
if strcmp(attrib,'LS') | strcmp(attrib,'CPCR') | strcmp(attrib,'CPCA') | strcmp(attrib,'MLR') | strcmp(attrib,'CDA') | strcmp(attrib,'CSIMPLS')|strcmp(attrib,'CSIMCA')
    oneplot=1; %only classical plots possible
elseif  (strcmp(attrib,'RAPCA') |strcmp(attrib,'MCDREG') | strcmp(attrib,'RPCR')|strcmp(attrib,'ROBPCA')| strcmp(attrib,'RSIMPLS')| strcmp(attrib,'MCDCOV')|strcmp(attrib,'LTS')|strcmp(attrib,'RDA')|strcmp(attrib,'RSIMCA'))  &  classicplots==0
    oneplot=2; %only robust plots possible
else 
    oneplot=0; %both robust and classical plot in the middle of the screen 
    bdwidth=5;
    topbdwidth=30;
    set(0,'Units','pixels');
    scnsize=get(0,'ScreenSize');
    pos1=[bdwidth, 1/3*scnsize(4)+bdwidth, scnsize(3)/2-2*bdwidth, scnsize(4)/2-(topbdwidth+bdwidth)];
    pos2=[pos1(1)+scnsize(3)/2, pos1(2), pos1(3), pos1(4)];
end

switch plotn
case 'all' %all possible plots
    close 
    switch attrib
        case {'MCDREG','MLR'}
            if oneplot==0
                figure('Position',pos1)
                regresdiagplot(Xdist,Rdist,cutoffX,cutoffR,k,multi,attrib,labsd,labresd)               
                figure('Position',pos2)
                regresdiagplot(Xdistcl,Rdistcl,cutoffXcl,cutoffRcl,k,multi,attribcl,labsd,labresd)               
            elseif oneplot==1
                figure
                regresdiagplot(Xdistcl,Rdistcl,cutoffXcl,cutoffRcl,k,multi,attribcl,labsd,labresd) %multi stond op 2
            elseif oneplot==2
                figure
                regresdiagplot(Xdist,Rdist,cutoffX,cutoffR,k,multi,attrib,labsd,labresd)
            end
        case {'CPCA', 'ROBPCA'}
            if oneplot==0
                figure('Position',pos1)
                screeplot(L,attrib)
                figure('Position',pos2)
                screeplot(Lcl,attribcl)
                pos1(2)=pos1(2)-40;
                pos2(2)=pos2(2)-40;
                figure('Position',pos1)
                scorediagplot(Xdist,OD,k,cutoffX,cutoffO,attrib,labsd,labod)  
                figure('Position',pos2)
                scorediagplot(Xdistcl,ODcl,k,cutoffXcl,cutoffOcl,attribcl,labsd,labod)
            elseif oneplot==1
                figure
                screeplot(Lcl,attrib)
                figure
                scorediagplot(Xdistcl,ODcl,k,cutoffXcl,cutoffOcl,attribcl,labsd,labod)          
            else
                figure
                screeplot(L,attrib)   
                figure
                scorediagplot(Xdist,OD,k,cutoffX,cutoffO,attrib,labsd,labod)
            end
        case {'RSIMPLS','CSIMPLS'}
                if oneplot==0
                    figure('Position',pos1)
                    scorediagplot(Xdist,OD,k,cutoffX,cutoffO,attrib,labsd,labod)   
                    figure('Position',pos2)
                    scorediagplot(Xdistcl,ODcl,k,cutoffXcl,cutoffOcl,attribcl,labsd,labod)
                    pos1(2)=pos1(2)-40;
                    pos2(2)=pos2(2)-40;
                    figure('Position',pos1)
                    regresdiagplot(Xdist,Rdist,cutoffX,cutoffR,k,multi,attribcl,labsd,labresd)   
                    figure('Position',pos2)
                    regresdiagplot(Xdistcl,Rdistcl,cutoffXcl,cutoffRcl,k,multi,attribcl,labsd,labresd)               
%3D plot is still highly memory consuming 08/04/04
%                     if exist('OD')
%                         Nihil=(OD <= 1.e-06);
%                         OD(Nihil)=0;
%                         help=OD;
%                     elseif exist('ODcl')
%                         Nihil=(ODcl <= 1.e-06);
%                         ODcl(Nihil)=0;
%                         helpcl=ODcl;
%                     end
                    
%                     pos1(2)=pos1(2)-80;
%                     pos2(2)=pos2(2)-80;
%                     if ~all(help) %in case k=rank(T)=r then OD = zero 
%                         figure('Position',pos1)
%                         regresdiagplot(Xdist,Rdist,cutoffX,cutoffR,k,multi,attrib,labsd,labresd)  
%                         figure('Position',pos2)
%                         regresdiagplot(Xdistcl,Rdistcl,cutoffXcl,cutoffRcl,k,multi,attribcl,labsd,labresd)
%                     else
%                         figure('Position',pos1)
%                         regresdiagplot3d(Xdist,OD,Rdist,cutoffX,cutoffO,cutoffR,k,attrib,multi,labsd,labod,labresd,0)
%                         title(['Robust 3D outlier map based on ',num2str(k),' LV'])
%                         figure('Position',pos2)
%                         regresdiagplot3d(Xdistcl,ODcl,Rdistcl,cutoffXcl,cutoffOcl,cutoffRcl,k,attribcl,multi,labsd,labod,labresd,0)
%                         title(['Classical 3D outlier map based on ',num2str(k),' LV'])
%                     end       
                elseif oneplot==1
                    figure
                    scorediagplot(Xdistcl,ODcl,k,cutoffXcl,cutoffOcl,attribcl,labsd,labod)         
                    figure
                    regresdiagplot(Xdistcl,Rdistcl,cutoffXcl,cutoffRcl,k,multi,attribcl,labsd,labresd)       
%3D plot is still highly memory consuming 08/04/04        
%                     if exist('ODcl')
%                         Nihil=(ODcl <= 1.e-06);
%                         ODcl(Nihil)=0;
%                         helpcl=ODcl;
%                     end                 
%                     figure
%                     if ~all(helpcl) %in case k=rank(T)=r then OD = zero 
%                         regresdiagplot(Rdistcl,Rdistcl,cutoffXcl,cutoffRcl,k,multi,'CPCR',labsd,labresd)
%                     else
%                         regresdiagplot3d(Xdistcl,ODcl,Rdistcl,cutoffXcl,cutoffOcl,cutoffRcl,k,attribcl,multi,labsd,labod,labresd,0)
%                         title(['Classical 3D outlier map based on ',num2str(k),' LV'])
%                     end       
                else
                    figure
                    scorediagplot(Xdist,OD,k,cutoffX,cutoffO,attrib,labsd,labod) 
                    figure
                    regresdiagplot(Xdist,Rdist,cutoffX,cutoffR,k,multi,attrib,labsd,labresd) 
%3D plot is still highly memory consuming 08/04/04
%                     if exist('OD')
%                         Nihil=(OD <= 1.e-06);
%                         OD(Nihil)=0;
%                         help=OD;
%                     end
%                     figure
%                     if ~all(help) %in case k=rank(T)=r then OD = zero % nog verfijnen is niet volledig nul is zeer klein (<10^-6)
%                         regresdiagplot(Xdist,Rdist,cutoffX,cutoffR,k,multi,attrib,labsd,labresd)  
%                     else
%                         regresdiagplot3d(Xdist,OD,Rdist,cutoffX,cutoffO,cutoffR,k,attrib,multi,labsd,labod,labresd,0)
%                         title(['Robust 3D outlier map based on ',num2str(k),' LV'])
%                     end       
                end
            case {'CPCR','RPCR'}
            if oneplot==0
                figure('Position',pos1)
                screeplot(L,attrib)
                figure('Position',pos2)
                screeplot(Lcl,attrib)
                pos1(2)=pos1(2)-40;
                pos2(2)=pos2(2)-40;
                figure('Position',pos1)
                scorediagplot(Xdist,OD,k,cutoffX,cutoffO,attrib,labsd,labod)   
                figure('Position',pos2)
                scorediagplot(Xdistcl,ODcl,k,cutoffXcl,cutoffOcl,attribcl,labsd,labod)
                pos1(2)=pos1(2)-80;
                pos2(2)=pos2(2)-80;
                figure('Position',pos1)
                regresdiagplot(Xdist,Rdist,cutoffX,cutoffR,k,multi,attribcl,labsd,labresd)   
                figure('Position',pos2)
                regresdiagplot(Xdistcl,Rdistcl,cutoffXcl,cutoffRcl,k,multi,attribcl,labsd,labresd)               
%3D plot is still highly memory consuming 08/04/04
%                 if exist('OD')
%                     Nihil=(OD <= 1.e-06);
%                     OD(Nihil)=0;
%                     help=OD;
%                 elseif exist('ODcl')
%                     Nihil=(ODcl <= 1.e-06);
%                     ODcl(Nihil)=0;
%                     helpcl=ODcl;
%                 end                       
%                 if ~all(help) %in case k=rank(T)=r then OD = zero 
%                     figure('Position',pos1)
%                     regresdiagplot(Xdist,Rdist,cutoffX,cutoffR,k,multi,attrib,labsd,labresd)  
%                     figure('Position',pos2)
%                     regresdiagplot(Xdistcl,Rdistcl,cutoffXcl,cutoffRcl,k,multi,attribcl,labsd,labresd)
%                 else
%                     figure('Position',pos1)
%                     regresdiagplot3d(Xdist,OD,Rdist,cutoffX,cutoffO,cutoffR,k,attrib,multi,labsd,labod,labresd,0)
%                     title(['Robust 3D outlier map based on ',num2str(k),' LV'])
%                     figure('Position',pos2)
%                     regresdiagplot3d(Xdistcl,ODcl,Rdistcl,cutoffXcl,cutoffOcl,cutoffRcl,k,attribcl,multi,labsd,labod,labresd,0)
%                     title(['Classical 3D outlier map based on ',num2str(k),' LV'])
%                 end       
            elseif oneplot==1
                figure
                screeplot(Lcl,attrib)    
                figure
                scorediagplot(Xdistcl,ODcl,k,cutoffXcl,cutoffOcl,attrib,labsd,labod)         
                figure
                regresdiagplot(Xdistcl,Rdistcl,cutoffXcl,cutoffRcl,k,multi,attribcl,labsd,labresd)               
%3D plot is still highly memory consuming 08/04/04
%                 if exist('ODcl')
%                     Nihil=(ODcl <= 1.e-06);
%                     ODcl(Nihil)=0;
%                     helpcl=ODcl;
%                 end                                
%                 figure
%                 if ~all(helpcl) %in case k=rank(T)=r then OD = zero 
%                     regresdiagplot(Xdistcl,Rdistcl,cutoffXcl,cutoffRcl,k,multi,'CPCR',labsd,labresd)
%                 else
%                     regresdiagplot3d(Xdistcl,ODcl,Rdistcl,cutoffXcl,cutoffOcl,cutoffRcl,k,attribcl,multi,labsd,labod,labresd,0)
%                     title(['Classical 3D outlier map based on ',num2str(k),' LV'])
%                 end       
            else
                figure
                screeplot(L,attrib)
                figure
                scorediagplot(Xdist,OD,k,cutoffX,cutoffO,attrib,labsd,labod) 
                figure
                regresdiagplot(Xdist,Rdist,cutoffX,cutoffR,k,multi,attrib,labsd,labresd) 
%3D plot is still highly memory consuming 08/04/04
%                 figure
%                 if exist('OD')
%                     Nihil=(OD <= 1.e-06);
%                     OD(Nihil)=0;
%                     help=OD;
%                 end                                 
%                 if ~all(help) %in case k=rank(T)=r then OD = zero 
%                     regresdiagplot(Xdist,Rdist,cutoffX,cutoffR,k,multi,attrib,labsd,labresd)  
%                 else
%                     regresdiagplot3d(Xdist,OD,Rdist,cutoffX,cutoffO,cutoffR,k,attrib,multi,labsd,labod,labresd,0)
%                     title(['Robust 3D outlier map based on ',num2str(k),' LV'])
%                 end       
            end
        case 'MCDCOV'
            if oneplot==0
                figure('Position',pos1)
                distplot(rd,cutoffR,attrib,labmcd)
                figure('Position',pos2)
                distplot(md,cutoffM,attribcl,labmcd)
                figure('Position',pos1)
                chiqqplot(rd,size(data,2),attrib)
                figure('Position',pos2)
                chiqqplot(md,size(data,2),attribcl)
                figure
                ddplot(md,rd,cutoffM,attrib,labmcd)
                figure 
                if size(data,2)~=2|data==0
                    axes
                    box on
                    text(0.05,0.5,'Tolerance ellipse is only available for bivariate data','color','r')
                else 
                    ellipsplot(center,covar,data,rd,labmcd)
                    set(findobj('color','r'),'color','m');
                    hold on
                    ellipsplot(centercl,covcl,data,md,labmcd)
                    set(findobj('color','r'),'color','b','linestyle',':');
                    hold off
                    legend_handles=[findobj('color','m'); findobj('color','b','linestyle',':')];
                    legend(legend_handles,'robust','classical','Location','best')
                end
            elseif oneplot==2
                figure
                distplot(rd,cutoffM,attrib,labmcd)
                figure
                chiqqplot(rd,size(data,2),attrib)
                figure
                ddplot(md,rd,cutoffM,attrib,labmcd)
                figure
                if size(data,2)~=2|data==0
                    axes
                    box on
                    text(0.05,0.5,'Tolerance ellipse is only available for bivariate data','color','r')
                else 
                    ellipsplot(center,covar,data,rd,labmcd,'MCDCOV')
                end
            end
        case {'LTS','LS'}
            if oneplot==0
                figure('Position',pos1)
                residualplot(fitted,standres,attrib,'Fitted values','Standardized LTS residual',lablts)
                figure('Position',pos2)
                residualplot(fittedcl,standrescl,attribcl,'Fitted value','Standardized LS residual',lablts) 
                pos1(2)=pos1(2)-40;
                pos2(2)=pos2(2)-40;                
                figure('Position',pos1)                
                residualplot(1:n,standres,attrib,'Index','Standardized LTS residual',lablts)
                figure('Position',pos2)
                residualplot(1:ncl,standrescl,attribcl,'Index','Standardized LS residual',lablts)       
                pos1(2)=pos1(2)-80;
                pos2(2)=pos2(2)-80;                 
                figure('Position',pos1)
                normqqplot(resid,attrib) 
                figure('Position',pos2)
                normqqplot(residcl,attribcl)         
                pos1(2)=pos1(2)-100;
                pos2(2)=pos2(2)-100;    
                figure('Position',pos1)
                if strcmp(attrib,'LTS')
                    regresdiagplot(Xdist,standres,cutoffX,cutoffY,k,0,attrib,labsd,labresd);
                else
                    regresdiagplot(Xdist,Rdist,cutoffX,cutoffR,k,1,attrib,labsd,labresd);
                end
                figure('Position',pos2)
                if strcmp(attribcl,'LS')
                    regresdiagplot(Xdistcl,standrescl,cutoffXcl,0,k,0,attribcl,labsd,labresd);
                else
                    regresdiagplot(Xdistcl,Rdistcl,cutoffXcl,cutoffRcl,k,1,attribcl,labsd,labresd);
                end
                pos1(2)=pos1(2)-120;
                pos2(2)=pos2(2)-120;  
                figure('Position',pos1)
                if size(x,2)~=1|x==0
                    axes
                    box on
                    text(0.05,0.5,'Scatter plot is only available for bivariate data','color','r')
                else 
                    lsscatter(x,y,fitted,attrib)
                end
                figure('Position',pos2)
                if size(xcl,2)~=1|xcl==0
                    axes
                    box on
                    text(0.05,0.5,'Scatter plot is only available for bivariate data','color','r')
                else
                    lsscatter(xcl,ycl,fittedcl,attribcl)
                end               
            elseif oneplot ==1
                figure
                residualplot(fittedcl,standrescl,attribcl,'Fitted value','Standardized LS residual',lablts) ;
                figure
                residualplot(1:ncl,standrescl,attribcl,'Index','Standardized LS residual',lablts) ;
                figure
                normqqplot(residcl,attrib);
                figure
                if strcmp(attrib,'LS')
                    regresdiagplot(Xdistcl,standrescl,cutoffXcl,0,k,0,attrib,labsd,labresd);
                else
                    regresdiagplot(Xdistcl,Rdistcl,cutoffXcl,cutoffRcl,k,1,attrib,labsd,labresd);
                end
                figure
                if size(xcl,2)~=1|xcl==0
                    axes
                    box on
                    text(0.05,0.5,'Scatter plot is only available for bivariate data','color','r')
                else
                    lsscatter(xcl,ycl,fittedcl,attribcl)
                end
              elseif oneplot==2
                figure
                residualplot(fitted,standres,attrib, 'Fitted value','Standardized LTS residual',lablts)
                figure
                residualplot(1:n,standres,attrib,'Index','Standardized LTS residual',lablts)
                figure
                normqqplot(resid,attrib)
                figure
                if strcmp(attrib,'LTS')
                    regresdiagplot(Xdist,standres,cutoffX,0,k,0,attrib,labsd,labresd);
                else
                    regresdiagplot(Xdist,Rdist,cutoffX,cutoffR,k,1,attrib,labsd,labresd);
                end
                figure
                if size(x,2)~=1|x==0
                    axes
                    box on
                    text(0.05,0.5,'Scatter plot is only available for bivariate data','color','r')
                else
                    lsscatter(x,y,fitted,attrib)
                end
            end
        case{'RDA','CDA'}
            if oneplot==0
                    figure('Position',pos1)
                    daplot(x,group,center,covar,class,method)
                    figure('Position',pos2)
                    daplot(xcl,groupcl,centercl,covcl,classcl,methodcl)
                elseif oneplot==1
                    figure
                    daplot(xcl,groupcl,centercl,covcl,classcl,methodcl)
                elseif oneplot==2
                    figure
                    daplot(x,group,center,covar,class,method)
                end
            case{'RSIMCA','CSIMCA'}
                if oneplot==0
                    figure('Position',pos1)
                    simcaplot(result)
                    figure('Position',pos2)
                    simcaplot(resultcl)
                elseif oneplot==1
                    figure
                    simcaplot(resultcl)
                elseif oneplot==2
                    figure
                    simcaplot(result)
                end    
        end
    choice=exitno;
case 'scree' %screeplot
    close 
    if oneplot==0
        figure('Position',pos1)
        screeplot(L,attrib)
        figure('Position',pos2)
        screeplot(Lcl,attribcl)
    elseif oneplot==1
        screeplot(Lcl,attribcl)
    else
        screeplot(L,attrib)
    end
case 'pcadiag' %diagnostic Scoreplot (Xdist-Odist)
    %close 
    if oneplot==0
        figure('Position',pos1)
        scorediagplot(Xdist,OD,k,cutoffX,cutoffO,attrib,labsd,labod)
        figure('Position',pos2)
        scorediagplot(Xdistcl,ODcl,k,cutoffXcl,cutoffOcl,attribcl,labsd,labod)
    elseif oneplot==1
        figure
        scorediagplot(Xdistcl,ODcl,k,cutoffXcl,cutoffOcl,attribcl,labsd,labod)
    else
        figure
        scorediagplot(Xdist,OD,k,cutoffX,cutoffO,attrib,labsd,labod)
    end
case '3ddiag' %3D-plot: highly memory consuming!!!!
    close 
    if oneplot==0
        if exist('OD')
            Nihil=(OD <= 1.e-06);
            OD(Nihil)=0;
            help=OD;
        elseif exist('ODcl')
            Nihil=(ODcl <= 1.e-06);
            ODcl(Nihil)=0;
            helpcl=ODcl;
        end
        if ~all(help)%in case k=rank(T)=r then OD = zero 
            figure('Position',pos1)
            regresdiagplot(Xdist,Rdist,cutoffX,cutoffR,k,multi,attrib,labsd,labresd)  
            figure('Position',pos2)
            regresdiagplot(Xdistcl,Rdistcl,cutoffXcl,cutoffRcl,k,multi,attribcl,labsd,labresd)
        else
            figure('Position',pos1)
            regresdiagplot3d(Xdist,OD,Rdist,cutoffX,cutoffO,cutoffR,k,attrib,multi,labsd,labod,labresd,0)
            title(['Robust 3D outlier map based on ',num2str(k),' LV'])
            figure('Position',pos2)
            regresdiagplot3d(Xdistcl,ODcl,Rdistcl,cutoffXcl,cutoffOcl,cutoffRcl,k,attribcl,multi,labsd,labod,labresd,0)
            title(['Classical 3D outlier map based on ',num2str(k),' LV'])
        end   
    elseif oneplot==1
        if exist('ODcl')
            Nihil=(ODcl <= 1.e-06);
            ODcl(Nihil)=0;
            helpcl=ODcl;
        end
        if ~all(helpcl) %in case k=rank(T)=r then OD = zero 
            regresdiagplot(Xdistcl,Rdistcl,cutoffXcl,cutoffRcl,k,multi,attribcl,labsd,labresd)
        else
            regresdiagplot3d(Xdistcl,ODcl,Rdistcl,cutoffXcl,cutoffOcl,cutoffRcl,k,attribcl,multi,labsd,labod,labresd,0)
            title(['Classical 3D outlier map based on ',num2str(k),' LV'])
        end       
    else
        if exist('OD')
            Nihil=(OD <= 1.e-06);
            OD(Nihil)=0;
            help=OD;
        end
        if ~all(help) %in case k=rank(T)=r then OD = zero 
            regresdiagplot(Xdist,Rdist,cutoffX,cutoffR,k,multi,attrib,labsd,labresd)  
        else
            regresdiagplot3d(Xdist,OD,Rdist,cutoffX,cutoffO,cutoffR,k,attrib,multi,labsd,labod,labresd,0)
            title(['Robust 3D outlier map based on ',num2str(k),' LV'])
        end  
    end
case {'Idist','robdist'} 
    close
    if oneplot==0
        figure('Position',pos1)
        distplot(rd,cutoffM,'MCDCOV',labmcd)
        figure('Position',pos2)
        distplot(md,cutoffM,'COV',labmcd)
        %elseif oneplot==1 not available since 'COV' is only called inside
        %the MCDCOV function
    elseif oneplot==2
        distplot(rd,cutoffM,'MCDCOV',labmcd)
    end
case 'qqmcd'
    close
    if oneplot==0
        figure('Position',pos1)
        chiqqplot(rd,size(data,2),'MCDCOV')
        figure('Position',pos2)
        chiqqplot(md,size(data,2),'COV')
    elseif oneplot==2
        chiqqplot(rd,size(data,2),'MCDCOV')
    end
case 'dd'
    close
    ddplot(md,rd,cutoffM,attrib,labmcd)
case 'ellipse'
    close
    if oneplot==0
        if size(data,2)~=2|data==0
            axes
            box on
            text(0.05,0.5,'Tolerance ellipse plot is only available for bivariate data','color','r')
        else
            ellipsplot(center,covar,data,rd,labmcd)
            set(findobj('color','r'),'color','m');
            hold on
            ellipsplot(centercl,covcl,data,md,labmcd)
            set(findobj('color','r'),'color','b','linestyle',':');
            hold off
            legend_handles=[findobj('color','m'); findobj('color','b','linestyle',':')];
            legend(legend_handles,'robust','classical','Location','best')
        end
    elseif oneplot==2
        if size(data,2)~=2|data==0
            axes
            box on
            text(0.05,0.5,'Tolerance ellipse plot is only available for bivariate data','color','r')
        else
            ellipsplot(center,covar,data,rd,labmcd)
        end
    end
case 'resfit'
    close
    if oneplot==0
        figure('Position',pos1)
        residualplot(fitted,standres,attrib,'Fitted value','Standardized LTS residual',lablts)
        figure('Position',pos2)
        residualplot(fittedcl,standrescl,attribcl,'Fitted value','Standardized LS residual',lablts) 
    elseif oneplot==1
        residualplot(fittedcl,standrescl,attribcl,'Fitted value','Standardized LS residual',lablts) 
    elseif oneplot==2
        residualplot(fitted,standres,attrib,'Fitted value','Standardized LTS residual',lablts)
    end
case 'resindex'
    close
    if oneplot==0
        figure('Position',pos1)
        residualplot(1:n,standres,attrib,'Index','Standardized LTS residual',lablts)
        figure('Position',pos2)
        residualplot(1:ncl,standrescl,attribcl,'Index','Standardized LS residual',lablts)  
    elseif oneplot==1
        residualplot(1:ncl,standrescl,attribcl,'Index','Standardized LS residual',lablts)  
    elseif oneplot==2
        residualplot(1:n,standres,attrib,'Index','Standardized LTS residual',lablts)      
    end
case 'qqlts'
    close
    if oneplot==0
        figure('Position',pos1)
        normqqplot(resid,'LTS')
        figure('Position',pos2)
        normqqplot(residcl,'LS')
    elseif oneplot==1
        normqqplot(residcl,'LS')
    elseif oneplot==2
        normqqplot(resid,'LTS')
    end
case 'regdiag'
    %close
    if oneplot==0
        figure('Position',pos1)
        if strcmp(attrib,'LTS')
            regresdiagplot(Xdist,standres,cutoffX,0,k,0,attrib,labsd,labresd);
        else
            regresdiagplot(Xdist,Rdist,cutoffX,cutoffR,k,multi,attrib,labsd,labresd);
        end
        figure('Position',pos2)
        if strcmp(attribcl,'LS')
            regresdiagplot(Xdistcl,standrescl,cutoffXcl,0,k,0,attribcl,labsd,labresd);
        else
            regresdiagplot(Xdistcl,Rdistcl,cutoffXcl,cutoffRcl,k,multi,attribcl,labsd,labresd);
        end
    elseif oneplot==1
        if strcmp(attrib,'LS')
            regresdiagplot(Xdistcl,standrescl,cutoffXcl,0,k,0,attrib,labsd,labresd);
        else
            regresdiagplot(Xdistcl,Rdistcl,cutoffXcl,cutoffRcl,k,multi,attrib,labsd,labresd);
        end
    elseif oneplot==2
        if strcmp(attrib,'LTS')
            regresdiagplot(Xdist,standres,cutoffX,0,k,0,attrib,labsd,labresd);
        else
            regresdiagplot(Xdist,Rdist,cutoffX,cutoffR,k,multi,attrib,labsd,labresd);
        end
    end
case 'scatter'
    close
    if oneplot==0
        figure('Position',pos1)
        if size(x,2)~=1|x==0
            axes
            box on
            text(0.05,0.5,'Scatter plot is only available for bivariate data','color','r')
        else 
            lsscatter(x,y,fitted,attrib)
        end
        figure('Position',pos2)
        if size(xcl,2)~=1|xcl==0
            axes
            box on
            text(0.05,0.5,'Scatter plot is only available for bivariate data','color','r')
        else
            lsscatter(xcl,ycl,fittedcl,attribcl)
        end
    elseif oneplot==1
        if size(xcl,2)~=1|xcl==0
            axes
            box on
            text(0.05,0.5,'Scatter plot is only available for bivariate data','color','r')
        else
            lsscatter(xcl,ycl,fittedcl,attribcl)
        end
    elseif oneplot==2
        if size(x,2)~=1|x==0
            axes
            box on
            text(0.05,0.5,'Scatter plot is only available for bivariate data','color','r')
        else
            lsscatter(x,y,fitted,attrib)
        end
    end
case 'da'
    close 
    if oneplot==0
        figure('Position',pos1)
        daplot(x,group,center,covar,class,method)
        figure('Position',pos2)
        daplot(xcl,groupcl,centercl,covcl,classcl,methodcl)
    elseif oneplot==1
        daplot(xcl,groupcl,centercl,covcl,classcl,methodcl)
    elseif oneplot==2
        daplot(x,group,center,covar,class,method)
    end 
case 'simca'
    close 
    if oneplot==0
        figure('Position',pos1)
        simcaplot(out)
        figure('Position',pos2)
        simcaplot(resultcl)
    elseif oneplot==1
        simcaplot(out)
    elseif oneplot==2
        simcaplot(out)
    end 
end

%---------------------------------------------------------------------------------------------------------------------------
function choice= menucov(out,attrib,exitno,labmcd,classic,classicplots)
choice=menu('Covariance plots: ','All','Index plot of the distances','Quantile-Quantile plot of the distances','Distance-distance plot',...
             'Tolerance ellipse (for bivariate data)','Exit');          
switch choice
    case 1
        plotn='all';
        whichplot(out,plotn,attrib,exitno,3,3,3,labmcd,3,classic,classicplots)
    case 2
        plotn='Idist';
        whichplot(out,'Idist',attrib,exitno,3,3,3,labmcd,3,classic,classicplots)
    case 3
        plotn='qqmcd';
        whichplot(out,plotn,attrib,exitno,3,3,3,labmcd,3,classic,classicplots)
    case 4
        plotn='dd';
        whichplot(out,plotn,attrib,exitno,3,3,3,labmcd,3,classic,classicplots)
    case 5
        plotn='ellipse';
        whichplot(out,plotn,attrib,exitno,3,3,3,labmcd,3,classic,classicplots)
    case 6
        choice=exitno;
    end
%---------------------------------------------------------------------------------------------------------------------------
function choice = menureg(out,attrib,exitno,labsd,labresd,classic,classicplots)
choice=menu('Regression Plots: ','All ','Regression outlier map', 'Exit ');
switch choice
    case 1
        plotn='all';
        whichplot(out,plotn,attrib,exitno,labsd,3,labresd,3,3,classic,classicplots)
    case 2
        plotn='regdiag';
        whichplot(out,plotn,attrib,exitno,labsd,3,labresd,3,3,classic,classicplots)
    case 3
        choice=exitno;
    end
%---------------------------------------------------------------------------------------------------------------------------
function choice =  menuscore(out,attrib,exitno,labsd,labod,classic,classicplots)
choice=menu('Score Plots: ','All ','Scree plot','Score outlier map', 'Exit ');
switch choice
    case 1
        plotn='all';
        whichplot(out,plotn,attrib,exitno,labsd,labod,3,3,3,classic,classicplots)
    case 2
        plotn='scree';
        whichplot(out,plotn,attrib,exitno,labsd,labod,3,3,3,classic,classicplots)
    case 3
        plotn='pcadiag';
        whichplot(out,plotn,attrib,exitno,labsd,labod,3,3,3,classic,classicplots)
    case 4
        choice=exitno;
    end
    
%---------------------------------------------------------------------------------------------------------------------------
function choice = menupls(out,attrib,exitno,labsd,labod,labresd,classic,classicplots)
choice=menu('PLS Plots: ','All ','Score outlier map',...
                     'Regression outlier map', '3D outlier map (Highly memory consuming)','Exit ');
switch choice
    case 1
        plotn='all';
        whichplot(out,plotn,attrib,exitno,labsd,labod,labresd,3,3,classic,classicplots)
    case 2
        plotn='pcadiag';
        whichplot(out,plotn,attrib,exitno,labsd,labod,labresd,3,3,classic,classicplots)  
    case 3
        plotn='regdiag';
        whichplot(out,plotn,attrib,exitno,labsd,labod,labresd,3,3,classic,classicplots) 
    case 4
        plotn='3ddiag';
        whichplot(out,plotn,attrib,exitno,labsd,labod,labresd,3,3,classic,classicplots)            
    case 5
        choice=exitno;
    end     
%---------------------------------------------------------------------------------------------------------------------------
function choice=menuscoreg(out,attrib,exitno,labsd,labod,labresd,classic,classicplots)
choice=menu('Score and Regression Plots: ','All ','Scree plot','Score outlier map',...
        'Regression outlier map', '3D outlier map  (Highly memory consuming)', 'Exit ');

switch choice
    case 1
        plotn='all';
        whichplot(out,plotn,attrib,exitno,labsd,labod,labresd,3,3,classic,classicplots)
    case 2
        plotn='scree';
        whichplot(out,plotn,attrib,exitno,labsd,labod,labresd,3,3,classic,classicplots)
    case 3
        plotn='pcadiag';
        whichplot(out,plotn,attrib,exitno,labsd,labod,labresd,3,3,classic,classicplots)              
    case 4
        plotn='regdiag';
        whichplot(out,plotn,attrib,exitno,labsd,labod,labresd,3,3,classic,classicplots)
    case 5
        plotn='3ddiag';
        whichplot(out,plotn,attrib,exitno,labsd,labod,labresd,3,3,classic,classicplots)
    case 6
        choice=exitno;
    end
%---------------------------------------------------------------------------------------------------------------------------
function choice= menuls(out,attrib,exitno,labsd,labresd,lablts,classic,classicplots)
choice=menu('Residual plots: ', 'All', 'Standardized residuals versus fitted values',...
             'Index plot of standardized residuals','Normal QQplot of residuals',...
             'Diagnostic plot of residuals versus robust distances',...
             'Scatter plot with regression line','Exit');         
switch choice
case 1
    plotn='all';
    whichplot(out,plotn,attrib,exitno,labsd,3,labresd,3,lablts,classic,classicplots)
case 2
    plotn='resfit';
    whichplot(out,plotn,attrib,exitno,3,3,3,3,lablts,classic,classicplots)
case 3
    plotn='resindex';
    whichplot(out,plotn,attrib,exitno,3,3,3,3,lablts,classic,classicplots)
case 4
    plotn='qqlts';
    whichplot(out,plotn,attrib,exitno,3,3,3,3,lablts,classic,classicplots)
case 5
    plotn='regdiag';
    whichplot(out,plotn,attrib,exitno,labsd,3,labresd,3,lablts,classic,classicplots)    
case 6
    plotn='scatter';
    whichplot(out,plotn,attrib,exitno,3,3,3,3,lablts,classic,classicplots)
case 7
    choice=exitno;
end   
%---------------------------------------------------------------------------------------------------------------------------
function choice = menuda(out,attrib,exitno,classic,classicplots)
choice=menu('Discriminant analysis: ','All ','Tolerance ellipses (for bivariate data)', 'Exit ');
switch choice
    case 1
        plotn='all';
        whichplot(out,plotn,attrib,exitno,3,3,3,3,3,classic,classicplots)
    case 2
        plotn='da';
        whichplot(out,plotn,attrib,exitno,3,3,3,3,3,classic,classicplots)
    case 3
        choice=exitno;
end
%------------------------------------------------------------------------------------------------------------------------------
function choice = menusimca(out,attrib,exitno,classic,classicplots)
choice=menu('SIMCA analysis: ','All ','Scatter plot', 'Exit ');
switch choice
    case 1
        plotn='all';
        whichplot(out,plotn,attrib,exitno,3,3,3,3,3,classic,classicplots)
    case 2
        plotn='simca';
        whichplot(out,plotn,attrib,exitno,3,3,3,3,3,classic,classicplots)
    case 3
        choice=exitno;
end
