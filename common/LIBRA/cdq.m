function result = cdq(x,y,c,varargin)

%CDQ computes Censored Depth Quantiles for regression, as described in
%
% Debruyne, M., Hubert, M., Portnoy, S., Vanden Branden, K. (2008),
% "Censored depth quantiles",
% Computational Statistics and Data Analysis, 52, 1604-1614.
%
% Required input arguments:
%    x : Data matrix of explanatory variables (also called 'regressors'). 
%        Rows of x represent observations, and columns represent variables.  
%        If the regression should contain an intercept, the x-matrix should
%        contain a column of ones.
%    y : A vector with n elements that contains the response variables.
%	 c : A vector that contains the indices of the censored observations.
%
% Optional input arguments:
%	 nrCan : Number of candidate hyperplanes to compute objective
%	     	 function at first grid point. The default is 500.
%	 nrCom : The objective function at each grid point after the first one is optimized over 
%            a set of nrCom hyperplanes having p-1 points in common with
%            the previously fitted hyperplane. The default is 100.
%            If nrCom = 0, the objective function at each grid point after the first one is optimized over 
%		     a fixed set of 'nrCan' hyperplanes.
%	 nrLam : Number of hyperplanes used in the computation of the regression depth of each hyperplane.
%            The default is 500.
%    grid  : Vector containing the gridpoints. The default is 0.05:0.05:0.95
%	 maxit : Maximum number of iterations at each grid point. The default
%	         is 4.
%
% I/O: result=cdq(x,y,c,'nrCan',nrCan,'nrCom',nrCom,'nrLam',nrLam,'grid',grid,'maxit',maxit);
%  The user should only give the input arguments that have to change their default value.
%  The name of the input arguments needs to be followed by their value.
%  The order of the input arguments is of no importance.
%
% The output of CDQ is a (g x p) matrix containing the regression quantiles,
% with g the number of grid points and p the number of regression parameters. 
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Michiel Debruyne
% Last Update: April 27, 2007.

%Handle defaults and optional user inputs.
nrCan=500;
nrCom=100;
nrLam=500;
grid=0.05:0.05:0.95;
maxit=4;
default=struct('nrCan',nrCan,'nrCom',nrCom,'nrLam',nrLam,'grid',grid,'maxit',maxit);
list=fieldnames(default);
options=default;
IN=length(list);
i=1;
%
if nargin==1
    error('At least 3 inputs are required, see help cdq.')
elseif nargin==2
    error('At least 3 inputs are required, see help cdq.')
end

if nargin>2
   
    nrCan=round(options.nrCan);
    nrLam=round(options.nrLam);
    nrCom=floor(options.nrCom);
    grid=sort(options.grid);
    maxit=round(options.maxit);

    for j=1:nargin-3 % searching the index of the accompanying field
        if rem(j,2)~=0 % fieldnames are placed on odd index
            if strcmp('nrCan',varargin{j})
                nrCan=varargin{j+1};
                if nrCan<1
                    error('nrCan should be positive, see help cdq.')
                end
            end
        end
        if rem(j,2)~=0 % fieldnames are placed on odd index
            if strcmp('nrCom',varargin{j})
                nrCom=varargin{j+1};
                if nrCom<1
                    error('nrCom should be positive, see help cdq.')
                end
            end
        end
        if rem(j,2)~=0 % fieldnames are placed on odd index
            if strcmp('nrLam',varargin{j})
                nrLam=varargin{j+1};
                if nrLam<1
                    error('nrLam should be positive, see help cdq.')
                end
            end
        end
        if rem(j,2)~=0 % fieldnames are placed on odd index
            if strcmp('grid',varargin{j})
                grid=varargin{j+1};
                if size(grid,1)>1 & size(grid,2)>1
                    error('grid should be a vector containing the grid points, not a matrix. See help cdq.')
                elseif size(grid,1)==1 & size(grid,2)==1
                    error('grid should be a vector containing the grid points, not a number. See help cdq.')
                elseif sum(grid>=1)+sum(grid<=0)
                    error('All entries of grid should be between 0 and 1. See help cdq.')
                end
            end
        end
        if rem(j,2)~=0 % fieldnames are placed on odd index
            if strcmp('maxit',varargin{j})
                maxit=varargin{j+1};
                if maxit<2
                    error('maxit should at least be 2, see help cdq.')
                end
            end
        end
    end
end

%Start main program
n=length(y);
p=size(x,2);
L=length(grid);

%Construct sets of hyperplanes.
M=nrCan+nrLam;
for i=1:M
   ok=0;
   while (ok==0)
      a=transpose(willcomb(n,p));
      ppu=x(a,:);   
      if (rank(ppu)==p)
         ok=1;
      end
   end
   beti=ppu\y(a);
   betakan(:,i)=beti;
end
setCan=betakan(:,1:nrCan);
setLam=betakan(:,nrCan+1 : nrCan+nrLam);

%tauh contains the tau-hats, defined for every crossed censored observation.
tauh=[];
crossedobs=[];

tau=grid(1);
beta(1,1)=tau;

%Determine the first estimate 
betaPrev=maxObj(x,y,tau,c,crossedobs,tauh,setCan,setLam);
tauPrev=tau;
beta(2:(p+1),1)=betaPrev;
l=2;
iterations=0;

while l<=L
   tau=grid(l);
   if nrCom~=0 & iterations==0
      setCom=bepKanU(x,y,betaPrev,nrCom);
      [betaNext,mm]=maxObj(x,y,tau,c,crossedobs,tauh,setCom,setLam);
   elseif nrCom==0
      betaNext=maxObj(x,y,tau,c,crossedobs,tauh,setCan,setLam);
   else
      [betaNext,mm]=maxObj(x,y,tau,c,crossedobs,tauh,setCom,setLam);
   end
   residuals=y-x*betaNext;
	
   
   %Good observations were crossed and still are.
   resn=residuals(crossedobs)<=0;
   goodObs=crossedobs(resn);
   
   %Bad observations were crossed but are not anymore.
   badObs=crossedobs(residuals(crossedobs)>0);
   
   %New ones are the ones that are crossed for the first time.
   notCross=comp(crossedobs,c);
   newObs=notCross(residuals(notCross)<=0);
   
 
   if iterations<maxit & ((length(badObs)>0)|(length(newObs)>0))
      
      %Retain the weights of the good ones.
      tauh=tauh(resn);
      crossedobs=goodObs;
      
      %Add the new ones.
      nrNew=length(newObs);
      crossedobs=[crossedobs newObs];
      
      %Their weight equals the previous grid point.
      tauh=[tauh repmat(tauPrev,[1,nrNew])];
      
      iterations=iterations+1;
      
   elseif iterations==maxit
       %No solution is found (zero in output).
        beta(:,end+1)=[tau; zeros(p,1)];
        iterations=0;
        residuals=y-x*betaNext;
        if ~isempty(tauh)
            tauh=tauh(tauh~=tauPrev);
            crossedobs=crossedobs(1:length(tauh));
        end
        l=l+1;
            
   %BetaNext is a good solution.
   else
       beta(:,end+1)=[tau;betaNext];
       iterations=0;
       betaPrev=betaNext;
       tauPrev=tau;
       l=l+1;
   end
end

result=beta(2:end,:)';

%function computing the maximum of the objective function

function[beta,maxobjf]=maxObj(x,y,tau,c,crossedobs,tauh,cand,setLam)

n=length(y);
nrCan=size(cand,2);
nrLam=size(setLam,2);

nc=comp(c,1:n);

i=1;
maxobjf=0;

while i<=nrCan
   
   betai=cand(:,i);
   ri=y-x*betai;
   
   minobjBetai=10^20;
   j=1;
   while j<=nrLam & minobjBetai>maxobjf
       
      %Take next direction.
      setLamj=betai-setLam(:,j);
      
      if ~isempty(find(abs(setLamj)>0.0000001))
         
         %Project.
         proj=x*setLamj;
         
         notCross=comp(crossedobs,c);
         Kcomp=[notCross nc];
         
         %Part objective function observations in Kcomp.
         riKcompos=ri(Kcomp)>=-0.0000001; projKcompos=proj(Kcomp)>0.0000001;
         riKcompne=ri(Kcomp)<=0.0000001; projKcompne=proj(Kcomp)<-0.0000001;
         obj1V=tau*sum(riKcompos.*projKcompne)+(1-tau)*sum(riKcompne.*projKcompos);
         obj1G=tau*sum(riKcompos.*projKcompos)+(1-tau)*sum(riKcompne.*projKcompne);
         
         if ~isempty(crossedobs)
            %Part objective function for pseudo-observations in crossed censured observations.
            xcl=x(crossedobs,:)*setLamj; ric=ri(crossedobs); xclp=xcl>0.0000001;
            xcln=xcl<-0.0000001; ricn=ric<=0.0000001;ricp=ric>=-0.0000001;
         	proposResnegV=tauh(xclp & ricn);
         	pronegResposV=tauh(xcln & ricp);
         	obj2V=tau*sum((tau-pronegResposV)./(1-pronegResposV))+(1-tau)*sum((tau-proposResnegV)./(1-proposResnegV));
            
            proposResnegG=tauh(xcln & ricn);
         	pronegResposG=tauh(xclp & ricp);
          	obj2G=tau*sum((tau-pronegResposG)./(1-pronegResposG))+(1-tau)*sum((tau-proposResnegG)./(1-proposResnegG));           

         	%Part objective function for pseudo-observations at infinity.
         	pronegV=tauh(xcl<-0.0000001);
            obj3V=tau*sum(1-(tau-pronegV)./(1-pronegV));
            
            pronegG=tauh(xcl>0.0000001);
         	obj3G=tau*sum(1-(tau-pronegG)./(1-pronegG));

      	else
            obj2V=0;
            obj2G=0;
            obj3V=0;
            obj3G=0;
        	end
         
         %Paste 3 parts.
         objG=obj1G+obj2G+obj3G;
         objV=obj1V+obj2V+obj3V;
                  
         minobjBetai=min([objG,objV,minobjBetai]);
      end
      j=j+1;
   end
   
   if minobjBetai>maxobjf
      beta=betai;
      maxobjf=minobjBetai;
   end
   i=i+1;
end

% ------------------------------------------------------------------------

function[thetann] = exchangeob(theta,Minv,zu,yu,zi,yi)

zit=transpose(zi);
u=Minv*zi;
w=-1/(1+zit*u);
thetan=theta-(yi-zit*theta)*w*u;
Mninv=Minv+w*u*transpose(u);

zut=transpose(zu);
un=Mninv*zu;
wn=-1/(1-zut*un);
thetann=thetan+(yu-zut*thetan)*wn*un;

% ------------------------------------------------------------------------

function [nr] = willcomb(n,p)

he=randperm(n);
nr=he(1:p);

% ------------------------------------------------------------------------

function[zComp]=comp(z,c)

n=length(z);
i=1;
while i<=n
   c=c(z(i)~=c);
   i=i+1;
end
zComp=c;

% ------------------------------------------------------------------------

function[setHyp]=bepKanU(x,y,betaNext,M)

betaNext;
n=size(x,1);
p=size(x,2);

residuals=y-x*betaNext;
a=abs(residuals)<=0.00001;

blq=cumsum(a);
if blq(end)>p
    a(blq>p)=0;
end

xa=x(a,:);
ya=y(a,1);
xb=x(a==0,:);
yb=y(a==0,1);
aa=size(xa);

setHyp=zeros(p,M);
setHyp(:,1)=betaNext;
minv=inv(transpose(xa)*xa);
for i=2:M
   j=ceil(rand(1,1)*p);
   k=ceil(rand(1,1)*(n-p));
   bb=size(xb);
   aa=size(xa);
   if rank([xa(1:(j-1),:);xa((j+1):end,:);xb(k,:)])==p
       gre=exchangeob(betaNext,minv,transpose(xa(j,:)),ya(j),transpose(xb(k,:)),yb(k,1));
       setHyp(:,i)=gre;
   else
       setHyp(:,i)=betaNext;
   end
end