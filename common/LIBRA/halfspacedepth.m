function result=halfspacedepth(u,v,x,y)

%HALFSPACEDEPTH computes the halfspace depth of a (two-dimensional) point theta 
% relative to a bivariate data set. 
%
% The algoritm is described in: 
%    Rousseeuw, P., Ruts, I. (1996),
%    "AS 307: Bivariate Location Depth",
%    Applied Statistics (JRSS-C), 45, 516-526.
% 
% Required input arguments:
%            u : first coordinate of the point theta
%            v : second coordinate of the point theta
%            x : vector containing the first coordinates of all the data
%                points
%            y : vector containing the second coordinates of all the data
%                points
%
% I/O: [result] = halfspacedepth(u,v,x,y);
%
% This function is part of the Matlab Library for Robust Analysis,
% available at: 
%           http://wis.kuleuven.be/stat/robust.html
%
% Written by Fabienne Verwerft
%
%
% Checking input
%
if not(length(x)==length(y))
    error('The vectors x and y must have the same length.')
end
if sum(isnan(x))>=1 || sum(isnan(y))>=1
    error('Missing values are not allowed')
end
%
% m is the number of data points that coincide with theta
%
eps=0.000001;
n=length(x);
norm=sqrt((x-u).^2 + (y-v).^2);
m=sum(norm<=eps);
ll=find(norm>eps);
norm=norm(ll);
xn=(x(ll)-u)./norm;
yn=(y(ll)-v)./norm;
%
% The vector containes the indices of the elements of x and y that satisfy
% the conditions.
%
k=find((abs(x(ll))>abs(y(ll))) & (xn>=0));
alfa(k)=asin(yn(k));
t=find(alfa(k)<0);
alfa(k(t))=2*pi+alfa(k(t));
k=find((abs(x(ll))>abs(y(ll))) & (xn<0));
alfa(k)=pi-asin(yn(k));
k=find((abs(x(ll))<=abs(y(ll))) & (yn>=0));
alfa(k)=acos(xn(k));
k=find((abs(x(ll))<=abs(y(ll))) & (yn<0));
alfa(k)=2*pi-acos(xn(k));
g=find(alfa>=(2*pi-eps));
alfa(g)=0;
%
nn=n-m;
if nn <= 1
    hdepth=m;
    result=hdepth;
    return
    % all data points coincide with theta
end
%
alfa=sort(alfa);
%
hoek=max(alfa(1)-alfa(nn)+2*pi,max(diff(alfa)));
if hoek > pi+eps
    hdepth=m;
    result=hdepth;
    return
    % hdepth=0 because theta lies outside the datacloud
end
%
% rotation around theta
% nu is the number of angles in the upper halfcircle
%
alfa=alfa-alfa(1);
nu=sum(alfa < (pi-eps));

if nu >= nn 
    hdepth=m;
    result=hdepth;
    return
    % hdepth=0 every angle in the upper halfcircle so theta outside the
    % datacloud.
end
%
% construction of the array F
%
beta=alfa+pi;
alfatwee=alfa+2*pi;
A=[alfa,alfatwee,beta];
Aindex(1:2*nn)=1;
Aindex((2*nn+1):3*nn)=0;
[As,Asin]=sort(A);
Aindexs=Aindex(Asin);
pp=cumsum(Aindexs);
juisten=find(Aindexs==0);
F=pp(juisten);
%
% Adjust the array F for the angles that coincide with beta.
%
gelijkab=intersect(find(diff(As)<=eps)+1,juisten);
betagelijka=Asin(gelijkab);
if length(gelijkab)>0
    for i=1:length(gelijkab)
        aantal=sum((As(Aindexs==1)+eps)<As(gelijkab(i)));
        F(betagelijka(i)-2*nn)=aantal;
    end
end
%
gindex=find(diff(alfa)==0)+1;
G=0:(nn-1);
%
% Adjust the array G for angles alfa who coincide
%
if length(gindex)>0
    for i=1:length(gindex)
        aantal=sum(alfa(1:gindex(i))<alfa(gindex(i)));
        G(gindex(i))=aantal;
    end
end
%
k=F-G;
numh=min(min(k,(nn-k)));
result=(numh+m);