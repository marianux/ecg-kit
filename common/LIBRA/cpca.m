function result=cpca(x,varargin)

%CPCA performs a classical principal components analysis 
% on the data set x. The loadings correspond with the eigenvectors of the
% classical covariance matrix of x.
% 
% Required input argument:
%      x : Data matrix (observations in the rows, variables in the
%          columns)
%
% Optional input argument: 
%      k : Number of principal components to compute. If k is missing, 
%          a screeplot is plotted allowing the selection of
%          the number of principal components.
%  plots : If equal to one (default), a menu is shown which allows to draw several plots,
%          such as a score outlier map. 
%          If 'plots' is equal to zero, all plots are suppressed.
%
% I/O: result=cpca(x,'k',2)
%
% The output is a structure containing 
%
%    result.P       : Classical loadings (eigenvectors of covariance matrix)
%    result.L       : Classical eigenvalues (eigenvalues of covariance matrix)
%    result.M       : Classical mean of the columns of x
%    result.T       : Classical scores 
%    result.k       : Number of (chosen) principal components 
%    result.sd      : Score (mahalanobis) distances within the classical PCA subspace
%    result.od      : Orthogonal distances to the classical PCA subspace 
%    result.cutoff  : Cutoff values for the score and orthogonal distances
%    result.flag     : The observations whose score distance is larger than result.cutoff.sd (==> result.flag.sd)
%                      or whose orthogonal distance is larger than result.cutoff.od (==> result.flag.od)
%                      can be considered as outliers and receive a flag equal to zero (result.flag.all).
%                      The regular observations receive a flag 1.
%    result.class   : 'CPCA' 
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Sabine Verboven 
% Last Update: 05/04/2003

default=struct('plots',1,'k',0);
list=fieldnames(default);
options=default;
IN=length(list);
i=1;   
counter=1;
if nargin>1
    %
    % placing inputfields in array of strings
    %
    for j=1:nargin-1
        if rem(j,2)~=0
            chklist{i}=varargin{j};
            i=i+1;
        end
    end 
    % Checking which default parameters have to be changed
    % and keep them in the structure 'options'.
    while counter<=IN 
        index=strmatch(list(counter,:),chklist,'exact');
        if ~isempty(index) % in case of similarity
            for j=1:nargin-2 % searching the index of the accompanying field
                if rem(j,2)~=0 % fieldnames are placed on odd index
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
end
[n,p]=size(x);
% using the optimal (least time consuming) algorithm
if n<p
	[P,T,L,r,Xc,out.M]=kernelEVD(x);
else
	[P,T,L,r,Xc,out.M]=classSVD(x);
end

if options.k==0
    screeplot(L,'CPCA');
    k=input('How many principal components would you like to retain? ');
else
    k=options.k;
end
out.k=k;
out.P=P(:,1:out.k);
out.T=T(:,1:out.k);
out.L=L(1:out.k)';

% Mahalanobis distance in classical PCA subspace
Tclas=Xc*out.P;
out.sd=sqrt(mahalanobis(Tclas,zeros(size(Tclas,2),1),'cov',out.L'));
out.cutoff.sd=sqrt(chi2inv(0.975,out.k));

% Orthogonal distances to classical PCA subspace
Xtilde=Tclas*out.P';
Cdiff=Xc-Xtilde;
for i=1:n
    out.od(i)=norm(Cdiff(i,:));
end
r=rank(x);
if k~=r
    m=mean(out.od.^(2/3));
    s=sqrt(var(out.od.^(2/3)));
    out.cutoff.od = sqrt(norminv(0.975,m,s)^3); 
else
    out.cutoff.od=0;
end

%Computing flags
if k~=p
    out.flag.od=(out.od<=out.cutoff.od)';
    out.flag.sd=(out.sd<=out.cutoff.sd)';
    out.flag.all=(out.flag.od)&(out.flag.sd);
else
    out.flag.od=(out.od<=out.cutoff.od)';
    out.flag.sd=(out.sd<=out.cutoff.sd)';
    out.flag.all=out.flag.sd;
end
out.class='CPCA';

result=struct('P',{out.P},'L',{out.L},'M',{out.M},'T',{out.T},'k',{out.k},...
    'sd', {out.sd'},'od',{out.od'},'cutoff',{out.cutoff},'flag',{out.flag},...
    'class',{out.class});

try
    if options.plots 
        makeplot(result)
    end
catch %output must be given even if plots are interrupted 
    %> delete(gcf) to get rid of the menu 
end