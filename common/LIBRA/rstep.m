function [S,P,t,kmax,med]= rstep(X,k,center,r); 

%RSTEP is an auxiliary function for 'rapca.m'.
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Created by Sabine Verboven and Mia Hubert (October 2000)
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
U=(Z -  repmat(median(Z,2),1,size(Z,2)))./(repmat(madc(Z')',1,size(Z,2)));
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

