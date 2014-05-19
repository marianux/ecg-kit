function result=rapca(data,varargin);

%RAPCA is a 'Reflection-based Algorithm for Principal Components Analysis'. 
% It is resistant to outliers in the data. The robust loadings are computed
% using projection-pursuit techniques and reflections. 
% Therefore RAPCA can be applied to both low and high-dimensional data sets.
% In low dimensions (at most 15), we recommend to use the MCD method instead
% (see mcdcov.m). 
%
% The RAPCA algorithm is described in 
%   Hubert, M., Rousseeuw, P.J., Verboven, S. (2002),
%   "A fast method for robust principal components with applications to chemometrics", 
%   Chemometrics and Intelligent Laboratory Systems, 60, 101-111.
%
% Required input arguments:
%         data : data matrix (observations in the rows, variables in the
%                columns)
%
% Optional input arguments:
%            k : number of principal components to compute
%    plots 0/1 : if equal to 1 a screeplot and an outlier map are drawn (default = 1)
%                else plots are suppressed 
%        labsd : the 'labsd' observations with largest score distance are
%                labeled on the outlier map (default = 3)
%        labod : the 'labod' observations with largest orthogonal distance are
%                labeled on the outlier map (default = 3)          
%	center 0/1 : if equal to 1 the data are centered around the L1-median (default = 1)
%                else the data are centered around the coordinatewise median
%                (not orthogonally equivariant, but faster)
%      classic : If equal to one, the classical PCA analysis will be performed
%                (see also cpca.m). (default = 0)
%
% If k is missing, or k = 0, a screeplot is drawn which allows you to select
% the number of principal components. If k = 0 and plots = 0, the algorithm itself 
% will determine the number of components. This is not recommended.
%
% I/O: result=rapca(x,'k',k,'plots',1,'labsd',3,'labod',3,'center',1,'classic',0);
%  The user should only give the input arguments that have to change their default value.
%  The name of the input arguments needs to be followed by their value.
%  The order of the input arguments is of no importance.
%  
% Examples: 
%    result=rapca(x,'k',3,'plots',0)
%    result=rapca(x,'labsd',5,'center',0)
%
% The output of RAPCA is a structure containing 
% 
%    result.P        : Robust loadings (eigenvectors)
%    result.L        : Robust eigenvalues       
%    result.M        : Robust center of the data
%    result.T        : Robust scores 
%    result.k        : Number of (chosen) principal components
%    result.sd       : Robust score distances within the robust PCA subspace
%    result.od       : Orthogonal distances to the robust PCA subspace 
%    result.cutoff   : Cutoff values for the robust score and orthogonal distances
%    result.flag     : The observations whose score distance is larger than result.cutoff.sd (==> result.flag.sd)
%                      or whose orthogonal distance is larger than result.cutoff.od (==> result.flag.od)
%                      can be considered as outliers and receive a flag equal to zero (result.flag.all).
%                      The regular observations receive a flag 1.
%    result.class    : 'RAPCA' 
%    result.classic  : If the input argument 'classic' is equal to one, this structure
%                      contains results of the classical PCA analysis (see also cpca.m). 
%
% Let n denote the number of observations, and p the number of original variables,
% then RAPCA finds a robust center (p x 1) of the data and a loading matrix P which 
% is (p x k) dimensional. Its columns are orthogonal and define a new coordinate
% system. The scores (n x k) are the coordinates of the centered observations with 
% respect to the loadings. The eigenvalues are the squared robust scales of the 
% observations projected on each of the loadings.
% Note that RAPCA also yields a robust covariance matrix (often singular) which
% can be computed as
%                         cov=result.P*result.L*result.P'
%
% The screeplot shows the eigenvalues and is helpful to select the number of 
% principal components.
% The outlier map visualizes the observations by plotting their orthogonal
% distance to the robust PCA subspace versus their robust distances 
% within the PCA subspace. This allows to classify the data points into 4 types:
% regular observations, good leverage points, bad leverage points and 
% orthogonal outliers. Remark that the RAPCA algorithm by construction passes
% through 'result.k' data points. The orthogonal distance of these data points is thus zero.
%
% The outlier map (or diagnostic plot) is described in
%    Hubert, M., Rousseeuw, P.J., Vanden Branden K. (2005),
%    "ROBPCA: a new approach to robust principal components analysis",  
%   Technometrics, 47, 64--79.
%  
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Sabine Verboven and Mia Hubert 
%
% Last Update: 23/12/2003


[n,p]=size(data);
counter=1;
default=struct('k',0,'center',1,'plots',1,'labsd',3,'labod',3,'h',floor(0.75*n),'classic',0);
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
end
k=options.k;
center=options.center;
plots=options.plots;

if k<0
    warning(['The number of principal components should be positive!']);
end

% First Step: classical SVD on the data 
% This step reduces the data space to the affine subspace 
% spanned by r=min(n-1,p) observations.
if n < p
   [loads,scores,lambda,r,centerX,clm]=kernelEVD(data);
else
   [loads,scores,lambda,r,centerX,clm]=classSVD(data);
end
X=scores;

% Second Step: Rstep on X 
% computes the robust eigenvectors and eigenvalues
[S,P,out.T,kmax,Rm]=rstep(X,k,center,r);
L=S'.^2;
out.P=loads*P;

if center==1
   % Retransforming the robust location to the original space
   out.M=clm+Rm*loads'; 
else
   out.M=median(data);
   datacentr=data-repmat(out.M,size(data,1),1);
   out.T=datacentr*out.P;
end

% Making screeplot to decide on the number of principal components
if plots==1 & k==0
    screeplot(L,'RAPCA');
    k=input(['How many principal components would you like to retain?\n']);
    k=max(min(k,kmax),1);
elseif plots==1    
    screeplot(L,'RAPCA');
    k=min(k,kmax);
elseif k~=0
    k=min(k,kmax)
else
    disp(['The number of principal components is defined by the algorithm.']); 
    disp(['It is set to ',num2str(kmax),'.']); 
    k=kmax;
end

% shrinking to k-dimensional subspace
out.P=out.P(:,1:k);
out.T=out.T(:,1:k);
out.L=L(1:k);
disp(['The outlier map is based on ',num2str(k),' principal component(s).'])
out.k=k; 
out.h=options.h;

% Computing distances 
% Robust score distances in robust PCA subspace
out.sd=sqrt(mahalanobis(out.T,zeros(size(out.T,2),1),'cov',out.L))';
out.cutoff.sd=sqrt(chi2inv(0.975,out.k));
% Orthogonal distances to robust PCA subspace
XRc=data-repmat(out.M,n,1);
Xtilde=out.T*out.P';
Rdiff=XRc-Xtilde;
for i=1:n
    out.od(i,1)=norm(Rdiff(i,:));
end
% Robust cutoff-value for the orthogonal distance
if k~=r
    [m,s]=unimcd(out.od.^(2/3),out.h);
    out.cutoff.od = sqrt(norminv(0.975,m,s).^3); 
else
    out.cutoff.od=0;
end
% Classical analysis
if options.classic==1
    out.classic.P=loads(:,1:out.k);
    out.classic.L=lambda(1:out.k);
    out.classic.M=clm;
    out.classic.T=scores(:,1:out.k);
    out.classic.k=out.k;
    % Mahalanobis distance in classical PCA subspace
    Tclas=centerX*loads(:,1:out.k);
    out.classic.sd=sqrt(mahalanobis(Tclas,zeros(size(Tclas,2),1),'cov',out.classic.L))';
    out.classic.cutoff.sd=sqrt(chi2inv(0.975,out.k));
    % Orthogonal distances to classical PCA subspace
    Xtilde=Tclas*loads(:,1:out.k)';
    Cdiff=centerX-Xtilde;
    for i=1:n
        out.classic.od(i,1)=norm(Cdiff(i,:));
    end
    % Classical cutoff-values
    if k~=r
        m=mean(out.classic.od.^(2/3));
        s=sqrt(var(out.classic.od.^(2/3)));
        out.classic.cutoff.od = sqrt(norminv(0.975,m,s)^3); 
    else
        out.classic.cutoff.od=0;
    end
    out.classic.cutoff.sd=sqrt(chi2inv(0.975,out.k));
    out.classic.flag.od=(out.classic.od<=out.classic.cutoff.od);
    out.classic.flag.sd=(out.classic.sd<=out.classic.cutoff.sd);
    out.classic.class='CPCA';
    out.classic.classic=1;
else
    out.classic=0;
end  

if k~=r
    out.flag.od=(out.od<=out.cutoff.od);
    out.flag.sd=(out.sd<=out.cutoff.sd);
    out.flag.all=(out.flag.od)&(out.flag.sd);
    if options.classic==1
        out.classic.flag.all=(out.classic.flag.od)&(out.classic.flag.sd);
    end
else
    out.flag.od=(out.od<=out.cutoff.od);
    out.flag.sd=(out.sd<=out.cutoff.sd);
    out.flag.all=out.flag.sd;
    if options.classic==1
        out.classic.flag.all=out.classic.flag.sd;
    end
end


% The output
result=struct('P',{out.P},'L',{out.L},'M',{out.M},'T',{out.T},'k',{out.k},...
    'sd', {out.sd},'od',{out.od},'cutoff',{out.cutoff},'flag',out.flag',...
    'class',{'RAPCA'},'classic',{out.classic});
% Making outlier map

try
    if plots & options.classic
        makeplot(result,'classic',1)
    elseif plots
        makeplot(result)
        %figure, scorediagplot(out.sd,out.od,out.k,out.cutoff.sd,out.cutoff.od,'RAPCA',options.labsd,options.labod)
    end
catch %output must be given even if plots are interrupted 
    %> delete(gcf) to get rid of the menu 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S,P,t,kmax,med]= rstep(X,k,center,r); 

%RSTEP: this is an auxiliary function for 'rapca.m'. 
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Created by Sabine Verboven and Mia Hubert 
% Part of the code is based on S-PLUS code from C. Croux.
%
% Last Update: 01/22/2002
% 

warning on;
if nargin<2
   k=0;
end
if nargin<3
   center=1;
end
if nargin<4
   r=rank(X);
end
[n,p]=size(X);
if k==0
   p1=min(floor(n/2),r);
else 
   p1=min([k,r,floor(n/2)]);
end

if k==0 | k > p1  
      disp(['The maximum number of principal components is ',num2str(p1),'.'])
      disp(['This is the minimum of (number of data points/2) and the rank of the data matrix.'])
end
S=zeros(p1,1);
V=zeros(p,p1);
switch center
case 0
   med=median(X);
   Xcentr=X-repmat(med,n,1); 
case 1
   med=l1median(X);
   Xcentr=X-repmat(med,n,1);
end
Xnewcentr=Xcentr;
kmax=0;
Transfo=eye(p);
for l=1:p1, 
   B=Xnewcentr;
   Bnorm=zeros(n,1);
   for i=1:n 
     Bnorm(i)=norm(B(i,:),2);
   end
   Bnormr=Bnorm(Bnorm > 1.e-12);
   B=B(Bnorm > 1.e-12,:);
   %Searching in directions A
   A=diag(1./Bnormr)*B;
   if size(Xnewcentr,2)==1 %case l=p1
      V(1:l-1,l)=0;
      V(l:p,l)=1;
      Vorigin(:,l)=Transfo*V(:,p1); 
      t=Xcentr*Vorigin(:,p1); %last step needs extraction of scale directly in p-dim space
      if n>40
         S(p1)=A_scale(t);
      else
         S(p1)=qnm(t);
      end
      kmax=kmax+1;
      break
   else
      Y=Xnewcentr*A'; %projected points in columns
   end
   if n>40
       s=A_scale(Y);
   else
       s=qnm(Y);
   end
   [c,vj]=sort(s); 
   j=vj(length(s)); 
   S(l)=s(j);
   if (S(1)/S(l) > 10^3) &(kmax<p1)
       l=p1+1;
       break
   else
       kmax=kmax+1;
   end
   if l==1 
      V(:,l)=A(j,:)';
   else
      V(1:l-1,l)=0;
      V(l:p,l)=A(j,:)';
   end
   % EigenVectors = columns of V
   % Constructing Transformation
   Base=eye(p-l+1);
   U=[];
   ndiff=norm(Base(:,1)-V(l:p,l),inf); %max norm of the normal vector
   if ndiff> 1.e-12
      if  V(l:p,l)'*Base(:,1) < 0
         V(l:p,l)=(-1)*V(l:p,l);
      end   
      u=(1./norm(Base(:,1)-V(l:p,l)))*(Base(:,1)-V(l:p,l));
      U=Base-2*repmat(u'*Base,p-l+1,1).*repmat(u,1,p-l+1);
   else
      U=Base; 
   end
   % Transforming eigenvectors to the original pxp dimensional space
   if l==1 
      Vorigin(:,l)=V(:,l);
      Transfo=U;
   else
      Edge=eye(p);
      Edge(l:p,l:p)=U;
      Vorigin(:,l)=Transfo*V(:,l);
      Transfo=Transfo*Edge;
   end
   Xnewcentr=Xnewcentr*U; 		%Reflection of data 
   Xnewcentr=removal(Xnewcentr,0,1);
end
[S,I]=greatsort(real(S(1:kmax)));
P=Vorigin(:,I);
t=Xcentr*P; 

%--------------------------------------------------------------------------
function [A_est]=A_scale(Z)

% A_SCALE calculates the A estimate of scale of the columns of Z
% 
% I/O: [A_est]=A_scale(Z); 
% 

Z=Z'; 
U=(Z - repmat(median(Z,2),1,size(Z,2)))./(repmat(madc(Z')',1,size(Z,2)));
[n,p]=size(U);
for i=1:n
   Ui=U(i,:);
   if any(isnan(Ui))
      scale(i)=0;
   else
      Zi=Z(i,:); 
      med=median(Zi); 
      m=madc(Zi-med);
      Zi=Zi(abs(Ui)<3.85);
      Ui=Ui(abs(Ui)<3.85);
      Ti=sqrt(sum((Ui.^2).*((3.85^2-Ui.^2).^4)))*sqrt(p)*0.9471*m;
      Ni=abs(sum((3.85^2-Ui.^2).*(3.85^2-5*(Ui.^2))));
      scale(i)=Ti/Ni;
   end
end
A_est=scale;
