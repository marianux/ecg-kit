function plotnumbers(x,y,number1,number2,ask,z,number3);

%PLOTNUMBERS marks the points (x,y,z) on a plot with their index number.
%
% Required input: 
%   x : x-coordinates of the plotted points
%   y : y-coordinates of the plotted points
%   
% Optional input:
%   number1 : number of  points with largest distance on x axes to identify
%   number2 : number of  points with largest distance on y axes to identify
%       ask :  0 = sorting by x
%              1 = sorting by y
%              2 = sorting by x and abs(y)
%              3 = sorting by x and z
%              4 = sorting by y and z
%              5 = sorting by x,y, and z.
%   z : z-coordinates of the plotted points
%   number3 :  number of points with largest distance on z axes to identify 
%
%I/O: plotnumbers(x,y,number1,number2,ask,z,number3);
%
%Written by S.Verboven
%Last revision: 22/12/2003

if nargin<2 
    error('Missing 2 input vectors x and y.')
end
if nargin==3
    number2=3;
    ask=0;
elseif nargin==2
    number1=3;
    number2=3;
    ask=0;
elseif nargin==4
    ask=0;
elseif nargin==6
    number3=3;    
end
n=length(x);
if ask==0
    % sort by x
    [ord,ind]=sort(x);
    ind=ind(n-number1+1:n)';
elseif ask==1
    % sort by y
    [ord,ind]=sort(y);
    ind=ind(n-number2+1:n)';
elseif ask==2
    % sort by x and y
    [ord1,ind1]=sort(x);
    [ord2,ind2]=sort(y);
    ind=zeros(number1+number2,1);
    ind(1:number1)=ind1(n-number1+1:n)';
    ind(number1+1:(number1+number2))=ind2(n-number2+1:n)';
elseif ask==3
    % sort by x and abs(y)
    [ord1,ind1]=sort(x);
    [ord2,ind2]=sort(abs(y));
    ind=zeros(number1+number2,1);
    ind(1:number1)=ind1(n-number1+1:n)';
    ind(number1+1:(number1+number2))=ind2(n-number2+1:n)';
elseif ask==4
    % sort by x and z
    [ord1,ind1]=sort(x);
    [ord2,ind2]=sort(z);
    ind=zeros(number1+number3,1);
    ind(1:number1)=ind1(n-number1+1:n)';
    ind(number1+1:(number1+number3))=ind2(n-number3+1:n)';
elseif ask==5
    % sort by y and z
    [ord1,ind1]=sort(y);
    [ord2,ind2]=sort(z);
    ind=zeros(number3+number2,1);
    ind(1:number2)=ind1(n-number2+1:n)';
    ind(number2+1:(number2+number3))=ind2(n-number3+1:n)';
elseif ask==6
    %sort by x, y and z
    [ord1,ind1]=sort(x);
    [ord2,ind2]=sort(y);
    [ord3,ind3]=sort(z);
    ind=zeros(number1+number2+number3,1);
    ind(1:number1)=ind1(n-number1+1:n)';
    ind(number1+1:number1+number2)=ind2(n-number2+1:n)';
    ind(number1+number2+1:number1+number2+number3)=ind3(n-number3+1:n)';
end
xrange=get(gca,'Xlim');
range=xrange(2)-xrange(1);
if nargin<6
    text(x(ind)+ range/50,y(ind),int2str(ind));
else
    text(x(ind)+ range/50,y(ind),z(ind),int2str(ind))
end
