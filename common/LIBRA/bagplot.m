function result=bagplot(x,varargin)

%BAGPLOT draws a bagplot, which is a generalisation of the univariate boxplot
% to bivariate data. The original bagplot is described in
%
% Rousseeuw, P.J., Ruts, I. and Tukey, J.W. (1999),
% "The bagplot: a bivariate boxplot", The American Statistician, 53, 382-387.
%
% The construction of this bagplot is based on the halfspacedepth (see also halfspacedepth.m).
% As the computation of the halfspacedepth is rather time-consuming, it is recommended at large datasets
% to perform the computations on a random subset of the data. The default size of the subset is 200, but this can be
% modified by the user.
%
% The bagplot can also be computed based on the adjusted outlyingness (see also adjustedoutlyingness.m).
% This method is described in:
%
% Hubert, M., and Van der Veeken, S. (2008),
%    "Outlier detection for skewed data", Journal of Chemometrics, to appear.
%
% For the up-to-date references, please consult the website:
%    http://wis.kuleuven.be/stat/robust.html
%
% The bagplot based on the adjusted outlyingness can be obtained by setting the optional input argument 'type' equal to 'ao'.
%
% Required input arguments:
%            x : bivariate data matrix (observations in the rows, variables in the
%                columns)
%
% Optional input arguments:
%   sizesubset : When drawing the bagplot based on the halfspacedepth, the size of the subset used to perform
%                the main computations.
%         type : To draw the bagplot based on the halfspacedepth, this parameter should be equal to 'hd' (default).
%                To draw the bagplot based on the adjusted outlyingness, it should be set to 'ao'.
%        plots : If equal to 1, a bagplot is drawn (default). If equal to zero, no plot is made.
%     colorbag : The color of the bag.
%   colorfence : The color of the fence.
%      databag : If this parameter is 1, the data within the bag are
%                plotted. If this parameter is set equal to 0, data points
%                are not displayed.
%    datafence : If this parameter is 1, the data within the fence are
%                plotted. If this parameter is set equal to 0, data points
%                are not displayed.
%
% I/O: result=bagplot(x,'sizesubset',200,'type','hd','colorbag',[],'colorfence',[],'databag',1,'datafence',1);
%  The user should only give the input arguments that have to change their default value.
%  The name of the input arguments needs to be followed by their value.
%  The order of the input arguments is of no importance.
%
% Examples:
%    result=bagplot(x,'colorfence',[0.2 0.2 0.5],'databag',0)
%    result=bagplot(x,'datafence',0,colorbag',[1 0 0])
%
%
% The output of bagplot is a structure containing
%
%    result.center       : center of the data. When 'type=hd', this corresponds with the Tukey median. When
%                          'type=ao', the point with smallest adjusted outlyingness 
%    result.type         : same of the input parameter type: 'hd' or 'ao'.
%    result.flag         : is 0 for outliers and equals 1 for regular points
%    result.datatype     : is 2 for observations in the bag, 1 for the observations in the
%                          fence and zero for outliers.
%    result.depth        : when 'type=hd' this is the halfspacedepth of each data point,
%                          when 'type=ao' the adjusted outlyingness.
%
% Written by Fabienne Verwerft on 25/05/2005
% Update and revision by Stephan Van der Veeken 10/12/2007
% Last revision: 11/02/2008, 25/02/2008

if nargin<1
    error('Input argument is undefined')
end

if size(x,1)<10
    error('At least 10 datapoints are needed for the calculation')
end

if size(x,2)~=2
    error('Data should be 2 dimensional')
end

S=x;
counter=1;
default=struct('colorbag',[0.6 0.6 1],'colorfence',[0.8 0.8 1],'sizesubset',200,...
    'databag',1,'datafence',1,'plots',1,'type','hd');
list=fieldnames(default);
result=default;
IN=length(list);
i=1;
%reading the user's input
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
    % and keep them in the structure 'result'.
    %
    while counter<=IN
        index=strmatch(list(counter,:),chklist,'exact');
        if ~isempty(index) %in case of similarity
            for j=1:nargin-1 %searching the index of the accompanying field
                if rem(j,2)~=0 %fieldnames are placed on odd index
                    if strcmp(chklist{index},varargin{j})
                        I=j;
                    end
                end
            end
            result=setfield(result,chklist{index},varargin{I+1});
            index=[];
        end
        counter=counter+1;
    end
end
colorbag=result.colorbag;
colorfence=result.colorfence;
databag=result.databag;
datafence=result.datafence;
plots=result.plots;
type=result.type;
sizesubset=result.sizesubset;

switch type
    case 'hd'
        s1=S(:,1);
        s2=S(:,2);
        Q=bagp(s1,s2,sizesubset);
        bag=Q.interpol;
        datatyp=Q.datatyp;
        tukm=Q.tukmed;
        depth=Q.depth;
    case 'ao'
        s1=S(:,1);
        s2=S(:,2);
        n=length(s1);
        Q = adjustedoutlyingness(S);
        D=[S,Q.outl,Q.flag,(1:n)'];
        P=sortrows(D,3);
        L=[P(:,1),P(:,2)];
        n=size(P,1);
        f=floor(n/2);
        g=sum(Q.flag);
        h=n-g;
        bag=L((1:f),:);
        hulp=[ones(f,1);2*ones(g-f,1);3*ones(h,1)];
        d1=[P,hulp];
        d2=sortrows(d1,5);
        datatyp=[d2(:,1),d2(:,2),d2(:,6)];
        tukm=L(1,:);
        depth=Q.outl;
end

i=find(datatyp(:,3)==1);
data1=datatyp(i,1:2);
i=find(datatyp(:,3)==2);
data2=datatyp(i,1:2);
data=[data1;data2];
i=find(datatyp(:,3)>=3);
outl=datatyp(i,1:2);
plak=[data;bag];
k=convhull(plak(:,1),plak(:,2));
whisk=plak(k,1:2);
if plots==1
    fill(whisk(:,1),whisk(:,2),colorfence,'LineStyle','none')
    hold on %This is essential in order to prevent that only the last executed plotting command is executed
    q=convhull(bag(:,1),bag(:,2));
    bagq=bag(q,1:2);
    fill(bagq(:,1),bagq(:,2),colorbag)%This line should be placed after the command plot(data...) to only %plot %the data outside the bag
    axis square
    if databag==1
        plot(data1(:,1),data1(:,2),'o','MarkerFaceColor','k','MarkerEdgeColor','k','Markersize',4)
    end
    if datafence==1
        plot(data2(:,1),data2(:,2),'o','MarkerFaceColor','k','MarkerEdgeColor','k','Markersize',4)
    end
    plot(outl(:,1),outl(:,2),'hk','MarkerFaceColor','k','Markersize',8)
    plot(tukm(1),tukm(2),'o','MarkerFaceColor','w','MarkerEdgeColor','w','MarkerSize',10);
    plot(tukm(1),tukm(2),'+k','Markersize',8)
    switch type
        case 'ao'
            title('bagplot based on adjusted outlyingness')
            set(gcf,'Name', 'Bagplot based on adjusted outlyingness', 'NumberTitle', 'off');
        case 'hd'
            title('bagplot based on halfspacedepth')
            set(gcf,'Name', 'Bagplot based on halfspacedepth', 'NumberTitle', 'off');
    end
    xmin=min(datatyp(:,1));
    xmax=max(datatyp(:,1));
    small=0.05*(xmax-xmin);
    xaxis1=xmin-small;
    xaxis2=xmax+small;
    ymin=min(datatyp(:,2));
    ymax=max(datatyp(:,2));
    small=0.05*(ymax-ymin);
    xaxis3=ymin-small;
    xaxis4=ymax+small;
    axis([xaxis1 xaxis2 xaxis3 xaxis4]);
    box on;
    hold off
end
Datatyp=datatyp(:,3);
flag=zeros(size(S,1),1);
for i=1:(size(S,1))
    if Datatyp(i)==3
        flag(i)=0;
    else
        flag(i)=1;
    end
end
datatype=zeros(size(S,1),1);
for i=1:(size(S,1))
    if Datatyp(i)==3
        datatype(i)=0;
    else if Datatyp(i)==2
            datatype(i)=1;
        else
            datatype(i)=2;
        end
    end
end
if type==1
    depth=zeros(size(S,1),1);
    for i=1:(size(S,1))
        depth(i,1)=halfspacedepth(s1(i,1),s2(i,1),s1,s2);
    end
end

result=struct('center',tukm,'depth',depth,'datatype',datatype,'flag',flag);

%-----------
function[pins]=rdraw(n,ntot)
pin=unidrnd(ntot,1,n);
pins=sort(pin);
uu=find(not(diff(pins)==0));
pins=pins(uu);
no=length(pins);
if no<n
    i=1;
    while i<=(n-no)
        in=unidrnd(ntot,1,1);
        if sum(pins==in)==0
            pins(no+i)=in;
            i=i+1;
        end
    end
end
pins=sort(pins);

%----
function [kount,ADK,empty]=isodepth(x,y,d,varargin)

% ISODEPTH is an algoritm that computes the depth region of a bivariate dataset
% corresponding to depth d.
% First, we have to check whether the data points are in general position. If not,
% a very small random number is added to each of the data points until the
% dataset comes in general position. All this is done in the m-file dithering.
% Then all special k-dividers must be found. The coordinates of the vertices
% of the depth region we are looking for are intersection points of these
% special k-dividers. So, consequently, every intersection point in turn has
% to be tested, for example by computing its depth (see halfspacedepth.m),
% to check whether it is a vertex of the depth region.
%
% The ISODEPTH algoritm is described in:
%    Ruts, I., Rousseeuw, P.J. (1996),
%    "Computing depth contours of bivariate point clouds",
%    Computational Statistics and Data Analysis, 23, 153-168.
%
% Required input arguments:
%            x : vector containing the first coordinates of all the data
%                points
%            y : vector containing the second coordinates of all the data
%                points
%            d : the depth of which the depth region has to be constructed
%
%
% I/O: [kount, ADK, empty]= isodepth(x,y,d);
%
% The output of isodepth is given by
%
%        kount : the total number of vertices of the depth region
%        ADK   : the coordinates of the vertices of the depth region
%        empty : logical value (1 if the depth region is empty, 0 if not)
%
% This function is part of the Matlab Library for Robust Analysis,
% available at:
%              http://wis.kuleuven.be/stat/robust.html
%
% Last Update: 29/04/2005

n=length(x);
eps=0.0000001;
%
% Checking input
%
if length(x)==1
    error('x is not a vector')
elseif not(length(x)==length(y))
    error('The vectors x and y must have the same length.')
end
if sum(isnan(x))>=1 || sum(isnan(y))>=1
    error('Missing values are not allowed')
end
if sum(x==x(1))==n
    error('All data points ly on a vertical line.')
elseif sum(y==y(1))==n
    error('All data points ly on a horizontal line.')
else
    R=corrcoef(x,y);
    if abs(R(1,2))==1
        error('All data points are collineair.')
    end
end
%
% Check whether the data is in general position. If not, add a very small random
% number to each of the data points.
%
[x,y, Index, angl, ind1,ind2]=dithering(x,y);

%
% Main part
%
if (n==1)&&(d==1)
    kount=n;
    ADK=[x,y];
    empty=0;
    return
end
%
if (d>floor(n/2))
    kount=0;
    ADK=0;
    empty=1;
    return
end
%
if n<=3
    kount=n;
    ADK=[x,y];
    empty=1;
    return
end
%
nrank(Index)=(1:n);
%
% Let the line rotate from zero to angle(1)
%
ncirq=Index;
kount=1;
halt=0;
M=length(angl);
if angl(1)>(pi/2)
    L=1;
    D1=ind1(L);
    IV1=nrank(D1);
    D2=ind2(L);
    IV2=nrank(D2);
    IV=ncirq(IV1);
    ncirq(IV1)=ncirq(IV2);
    ncirq(IV2)=IV;
    IV=IV1;
    nrank(D1)=IV2;
    nrank(D2)=IV;
    %
    if ((IV1==d) && (IV2==(d+1)))||((IV2==d) && (IV1==(d+1)))||((IV1==(n-d)) && (IV2==(n-d+1)))||((IV2==(n-d)) && (IV1==(n-d+1)))
        if angl(L)<(pi/2)
            dum=angl(L)+(pi/2);
        else
            dum=angl(L)-(pi/2);
        end
        if (IV1==d && IV2==(d+1))||(IV2==d && IV1==(d+1))
            if dum<=(pi/2)
                alfa(kount)=angl(L)+pi;
            else
                alfa(kount)=angl(L);
            end
        end
        if or((IV1==(n-d) && IV2==(n-d+1)),(IV2==(n-d) && IV1==(n-d+1)))
            if dum<=(pi/2)
                alfa(kount)=angl(L);
            else
                alfa(kount)=angl(L)+pi;
            end
        end
        kand1(kount)=ind1(L);
        kand2(kount)=ind2(L);
        D(kount)=sin(alfa(kount))*x(ind1(L))-cos(alfa(kount))*y(ind1(L));
        kount=kount+1;
    end
    halt=1;
end
%
L=2;
stay=1;
% jflag keeps track of which angle we have to test next
while stay==1
    stay=0;
    kontrol=0;
    if (pi<=(angl(L)+(pi/2))) && ((angl(L)-(pi/2))< angl(1))
        D1=ind1(L);
        IV1=nrank(D1);
        D2=ind2(L);
        IV2=nrank(D2);
        IV=ncirq(IV1);
        ncirq(IV1)=ncirq(IV2);
        ncirq(IV2)=IV;
        IV=IV1;
        nrank(D1)=IV2;
        nrank(D2)=IV;
        %
        if ((IV1==d) && (IV2==(d+1)))||((IV2==d) && (IV1==(d+1)))||((IV1==(n-d)) && (IV2==(n-d+1)))||((IV2==(n-d)) && (IV1==(n-d+1)))
            if angl(L)<(pi/2)
                dum=angl(L)+(pi/2);
            else
                dum=angl(L)-(pi/2);
            end
            if (IV1==d && IV2==(d+1))||(IV2==d && IV1==(d+1))
                if dum<=(pi/2)
                    alfa(kount)=angl(L)+pi;
                else
                    alfa(kount)=angl(L);
                end
            end
            if or((IV1==(n-d) && IV2==(n-d+1)),(IV2==(n-d) && IV1==(n-d+1)))
                if dum<=(pi/2)
                    alfa(kount)=angl(L);
                else
                    alfa(kount)=angl(L)+pi;
                end
            end
            kand1(kount)=ind1(L);
            kand2(kount)=ind2(L);
            D(kount)=sin(alfa(kount))*x(ind1(L))-cos(alfa(kount))*y(ind1(L));
            kount=kount+1;
        end
        kontrol=1;
    end
    L=L+1;
    if kontrol==1
        halt=1;
    end
    if (L==(M+1)) && (kontrol==1)
        jflag=1;
        stay=2;
    end
    if not(stay==2)
        if ((halt==1)&&(kontrol==0))||(L==(M+1))
            stay=3;
        else
            stay=1;
        end
    end
end
if not(stay==2)
    if (L>1)
        jflag=L-1;
    else
        jflag=M;
    end
end
%
halt2=0;
if not(stay==2)
    J=0;
    %
    % If the first switch didnt occur between 0 and the angle angl(1) look for it
    % between the following angles.
    %
    stay2=1;
    if (L==M+1) && (kontrol==0)
        halt=0;
        halt2=0;
        J=J+1;
        if J==(M+1)
            J=1;
        end
        L=J+1;
        if L==(M+1)
            L=1;
        end
        while stay2==1
            stay2=0;
            kontrol=0;
            if (angl(L)+pi/2)<pi
                ang1=angl(L)+pi/2;
            else
                ang1=angl(L)-pi/2;
            end
            if J==M
                jj=1;
                if halt2==0
                    angl(1)=angl(1)+pi;
                end
            else
                jj=J+1;
            end
            if (angl(J)<=ang1) && (ang1<angl(jj))
                if angl(1)>pi
                    angl(1)=angl(1)-pi;
                end
                D1=ind1(L);
                IV1=nrank(D1);
                D2=ind2(L);
                IV2=nrank(D2);
                IV=ncirq(IV1);
                ncirq(IV1)=ncirq(IV2);
                ncirq(IV2)=IV;
                IV=IV1;
                nrank(D1)=IV2;
                nrank(D2)=IV;
                %
                if ((IV1==d) && (IV2==(d+1)))||((IV2==d) && (IV1==(d+1)))||((IV1==(n-d)) && (IV2==(n-d+1)))||((IV2==(n-d)) && (IV1==(n-d+1)))
                    if angl(L)<(pi/2)
                        dum=angl(L)+(pi/2);
                    else
                        dum=angl(L)-(pi/2);
                    end
                    if (IV1==d && IV2==(d+1))||(IV2==d && IV1==(d+1))
                        if dum<=(pi/2)
                            alfa(kount)=angl(L)+pi;
                        else
                            alfa(kount)=angl(L);
                        end
                    end
                    if or((IV1==(n-d) && IV2==(n-d+1)),(IV2==(n-d) && IV1==(n-d+1)))
                        if dum<=(pi/2)
                            alfa(kount)=angl(L);
                        else
                            alfa(kount)=angl(L)+pi;
                        end
                    end
                    kand1(kount)=ind1(L);
                    kand2(kount)=ind2(L);
                    D(kount)=sin(alfa(kount))*x(ind1(L))-cos(alfa(kount))*y(ind1(L));
                    kount=kount+1;
                end
                kontrol=1;
            end
            if angl(1)>pi
                angl(1)=angl(1)-pi;
            end
            if L==M
                L=1;
            else
                L=L+1;
            end
            if kontrol==1
                halt=1;
            end
            if (halt==1)&&(kontrol==0)
                if halt2==1
                    stay2=2;
                end
                if not(stay2==2)
                    if L>1
                        jflag=L-1;
                    else
                        jflag=M;
                    end
                    stay2=0;
                end
            else
                if L==jj
                    if jj==1
                        halt2=1;
                    end
                    J=J+1;
                    if J==(M+1)
                        J=1;
                    end
                    L=J+1;
                    if L==(M+1)
                        L=1;
                    end
                    stay2=1;
                else
                    stay2=1;
                end
            end
        end
    end
end
%
if not(stay2==2)
    %
    % The first switch has occurred. Now start looking for the next ones,
    % between the following angles.
    %
    for i=(J+1):(M-1)
        L=jflag;
        stay=1;
        while stay==1
            stay=0;
            kontrol=0;
            if ((angl(L)+pi/2)<pi)
                ang1=angl(L)+pi/2;
            else
                ang1=angl(L)-pi/2;
            end
            if (angl(i)<=ang1)&&(ang1<angl(i+1))
                D1=ind1(L);
                IV1=nrank(D1);
                D2=ind2(L);
                IV2=nrank(D2);
                IV=ncirq(IV1);
                ncirq(IV1)=ncirq(IV2);
                ncirq(IV2)=IV;
                IV=IV1;
                nrank(D1)=IV2;
                nrank(D2)=IV;
                %
                if ((IV1==d) && (IV2==(d+1)))||((IV2==d) && (IV1==(d+1)))||((IV1==(n-d)) && (IV2==(n-d+1)))||((IV2==(n-d)) && (IV1==(n-d+1)))
                    if angl(L)<(pi/2)
                        dum=angl(L)+(pi/2);
                    else
                        dum=angl(L)-(pi/2);
                    end
                    if (IV1==d && IV2==(d+1))||(IV2==d && IV1==(d+1))
                        if dum<=(pi/2)
                            alfa(kount)=angl(L)+pi;
                        else
                            alfa(kount)=angl(L);
                        end
                    end
                    if or((IV1==(n-d) && IV2==(n-d+1)),(IV2==(n-d) && IV1==(n-d+1)))
                        if dum<=(pi/2)
                            alfa(kount)=angl(L);
                        else
                            alfa(kount)=angl(L)+pi;
                        end
                    end
                    kand1(kount)=ind1(L);
                    kand2(kount)=ind2(L);
                    D(kount)=sin(alfa(kount))*x(ind1(L))-cos(alfa(kount))*y(ind1(L));
                    kount=kount+1;
                end
                kontrol=1;
            end
            if kontrol==0
                jflag=L;
            else
                if not(L==M)
                    L=L+1;
                else
                    L=1;
                end
                stay=1;
            end
        end
    end
    L=jflag;
    %
    % Finally, look for necessary switches between the last angle and zero.
    %
    stay=1;
    while stay==1
        kontrol=0;
        stay=0;
        if (angl(L)+pi/2)<pi
            ang1=angl(L)+pi/2;
        else
            ang1=angl(L)-pi/2;
        end
        if (angl(M)<=ang1)&&(ang1<pi)
            D1=ind1(L);
            IV1=nrank(D1);
            D2=ind2(L);
            IV2=nrank(D2);
            IV=ncirq(IV1);
            ncirq(IV1)=ncirq(IV2);
            ncirq(IV2)=IV;
            IV=IV1;
            nrank(D1)=IV2;
            nrank(D2)=IV;
            %
            if ((IV1==d) && (IV2==(d+1)))||((IV2==d) && (IV1==(d+1)))||((IV1==(n-d)) && (IV2==(n-d+1)))||((IV2==(n-d)) && (IV1==(n-d+1)))
                if angl(L)<(pi/2)
                    dum=angl(L)+(pi/2);
                else
                    dum=angl(L)-(pi/2);
                end
                if (IV1==d && IV2==(d+1))||(IV2==d && IV1==(d+1))
                    if dum<=(pi/2)
                        alfa(kount)=angl(L)+pi;
                    else
                        alfa(kount)=angl(L);
                    end
                end
                if or((IV1==(n-d) && IV2==(n-d+1)),(IV2==(n-d) && IV1==(n-d+1)))
                    if dum<=(pi/2)
                        alfa(kount)=angl(L);
                    else
                        alfa(kount)=angl(L)+pi;
                    end
                end
                kand1(kount)=ind1(L);
                kand2(kount)=ind2(L);
                D(kount)=sin(alfa(kount))*x(ind1(L))-cos(alfa(kount))*y(ind1(L));
                kount=kount+1;
            end
            kontrol=1;
        end
        if kontrol==1
            if not(L==M)
                L=L+1;
            else
                L=1;
            end
            stay=1;
        end
    end
end
num=kount-1; % num is the total number of special k-dividers
%
% Sort the num special k-dividers. Permute kand1, kand2 and D in the same
% way.
%
[alfa,In]=sort(alfa);
kand1=kand1(In);
kand2=kand2(In);
D=D(In);
%
IW1=1;
IW2=2;
Jfull=0;
NDK=0;
stay2=1;
while stay2==1
    stay2=0;
    ndata=0;
    %
    % Compute the intersection point.
    %
    while abs(-sin(alfa(IW2))*cos(alfa(IW1))+sin(alfa(IW1))*cos(alfa(IW2)))<eps
        IW2=IW2+1;
        ndata=0;
        if IW2==(num+1)
            IW2=1;
        end
    end
    %
    xcord=(cos(alfa(IW2))*D(IW1)-cos(alfa(IW1))*D(IW2))/(-sin(alfa(IW2))*cos(alfa(IW1))+sin(alfa(IW1))*cos(alfa(IW2)));
    ycord=(-sin(alfa(IW2))*D(IW1)+sin(alfa(IW1))*D(IW2))/(-sin(alfa(IW1))*cos(alfa(IW2))+sin(alfa(IW2))*cos(alfa(IW1)));
    %
    % Test whether the intersection point is a data point. If so,
    % adjust IW1 and IW2.
    %
    if or(kand1(IW1)==kand1(IW2),kand1(IW1)==kand2(IW2))
        ndata=kand1(IW1);
    end
    if or(kand2(IW1)==kand1(IW2),kand2(IW1)==kand2(IW2))
        ndata=kand2(IW1);
    end
    if not(ndata==0)
        iv=0;
        stay=1;
        while stay==1
            stay=0;
            next=IW2+1;
            iv=iv+1;
            if next==(num+1)
                next=1;
            end
            if not(next==IW1)
                if or(ndata==kand1(next),ndata==kand2(next))
                    IW2=IW2+1;
                    if (IW2==(num+1))
                        IW2=1;
                    end
                    stay=1;
                end
            end
        end
        if iv==(num-1)
            kount=1;
            ADK=[x(ndata),y(ndata)];
            empty=0;
            return
        end
    end
    if IW2==num
        kon=1;
    else
        kon=IW2+1;
    end
    if kon==IW1
        kon=kon+1;
    end
    if kon==(num+1)
        kon=1;
    end
    %
    % Test whether the intersection point lies to the left of the special
    % k-divider which corresponds to alfa(kon). If so, compute its depth.
    %
    stay3=1;
    stay4=1;
    if (sin(alfa(kon))*xcord-cos(alfa(kon))*ycord-D(kon))<=eps
        hdep1=halfspacedepth(xcord,ycord,x,y);
        if hdep1==d
            NDK=1;
        else
            hdep2=halfspacedepth(xcord-0.000001,ycord-0.000001,x,y);
            hdep3=halfspacedepth(xcord+0.000001,ycord+0.000001,x,y);
            hdep4=halfspacedepth(xcord-0.000001,ycord+0.000001,x,y);
            hdep5=halfspacedepth(xcord+0.000001,ycord-0.000001,x,y);
            hdepvector=[hdep1;hdep2;hdep3;hdep4;hdep5];
            if (NDK==0)&&(sum(hdepvector>=d)>=1)
                NDK=1;
            end
            if (hdep1<d)&&(hdep2<d)&&(hdep3<d)&&(hdep4<d)&&(hdep5<d)&&(NDK==1)
                %
                % The intersection point is not the correct one, try the next
                % special k-divider.
                %
                IW2=IW2+1;
                if IW2==(num+1)
                    IW2=1;
                end
                stay2=1;
            end
        end
        if not(stay2==1)
            %
            % Store IW2 and IW1 in kornr. If kornr has already been filled,
            % check wether we have encountered this intersection point before.
            %
            if (IW2>IW1)&&(Jfull==0)
                kornr(IW1:(IW2-1),1)=kand1(IW1);
                kornr(IW1:(IW2-1),2)=kand2(IW1);
                kornr(IW1:(IW2-1),3)=kand1(IW2);
                kornr(IW1:(IW2-1),4)=kand2(IW2);
            else
                if IW2>IW1
                    i=IW1;
                    stay3=1;
                    while stay3==1
                        if (kornr(i,1)==kand1(IW1))&&(kornr(i,2)==kand2(IW1))&&(kornr(i,3)==kand1(IW2))&&(kornr(i,4)==kand2(IW2))
                            stay3=0;
                        else
                            m1=(y(kornr(i,2))-y(kornr(i,1)))/(x(kornr(i,2))-x(kornr(i,1)));
                            m2=(y(kornr(i,4))-y(kornr(i,3)))/(x(kornr(i,4))-x(kornr(i,3)));
                            if not(m1==m2)
                                xcord1=(m1*x(kornr(i,1))-y(kornr(i,1))-m2*x(kornr(i,3))-y(kornr(i,3)))/(m1-m2);
                                ycord1=(m2*(m1*x(kornr(i,1))-y(kornr(i,1)))-m1*(m2*x(kornr(i,3))-y(kornr(i,3))))/(m1-m2);
                            end
                            if (abs(xcord1-xcord)<eps)&&(abs(ycord1-ycord)<eps)
                                stay3=0;
                            end
                            if stay3==1
                                kornr(i,1)=kand1(IW1);
                                kornr(i,2)=kand2(IW1);
                                kornr(i,3)=kand1(IW2);
                                kornr(i,4)=kand2(IW2);
                            end
                        end
                        if stay3==1
                            i=i+1;
                            if i==IW2
                                stay3=2;
                                i=i-1;
                            end
                        end
                    end
                else
                    Jfull=1;
                    kornr(IW1:num,1)=kand1(IW1);
                    kornr(IW1:num,2)=kand2(IW1);
                    kornr(IW1:num,3)=kand1(IW2);
                    kornr(IW1:num,4)=kand2(IW2);
                    i=1;
                    stay4=1;
                    if IW2==1
                        stay4=3;
                    end
                    while stay4==1
                        if (kornr(i,1)==kand1(IW1))&&(kornr(i,2)==kand2(IW1))&&(kornr(i,3)==kand1(IW2))&&(kornr(i,4)==kand2(IW2))
                            stay4=0;
                        else
                            m1=(y(kornr(i,2))-y(kornr(i,1)))/(x(kornr(i,2))-x(kornr(i,1)));
                            warning off MATLAB:divideByZero
                            m2=(y(kornr(i,4))-y(kornr(i,3)))/(x(kornr(i,4))-x(kornr(i,3)));
                            warning off MATLAB:divideByZero
                            if not(m1==m2)
                                xcord1=(m1*x(kornr(i,1))-y(kornr(i,1))-m2*x(kornr(i,3))-y(kornr(i,3)))/(m1-m2);
                                ycord1=(m2*(m1*x(kornr(i,1))-y(kornr(i,1)))-m1*(m2*x(kornr(i,3))-y(kornr(i,3))))/(m1-m2);
                            end
                            if (abs(xcord1-xcord)<=eps)&&(abs(ycord1-ycord)<=eps)
                                stay4=0;
                            end
                            if stay4==1
                                kornr(i,1)=kand1(IW1);
                                kornr(i,2)=kand2(IW1);
                                kornr(i,3)=kand1(IW2);
                                kornr(i,4)=kand2(IW2);
                            end
                        end
                        if stay4==1
                            i=i+1;
                            if i==IW2
                                i=i-1;
                                stay4=2;
                            end
                        end
                    end
                end
            end
        end
    elseif (stay3>0)&&(stay4>0)&&not(stay2==1)
        %
        % The intersection point is not the correct one, try the next
        % special k-divider.
        %
        IW2=IW2+1;
        if IW2==(num+1)
            IW2=1;
        end
        stay2=1;
    end
    %
    % Look for the next vertex of the convex figure.
    %
    if (stay3>0)&&(stay4>0)&&not(stay2==1)
        IW1=IW2;
        IW2=IW2+1;
        if IW2==(num+1)
            IW2=1;
        end
        stay2=1;
    end
end
%
% Scan kornr and ascribe the coordinates of the vertices to the variable
% ADK.
%
kount=0;
%
if NDK==0
    %
    % The requested depth region is empty
    %
    ADK=0;
    empty=1;
    return
else
    empty=0;
end
%
i=1;
E1=y(kornr(i,2))-y(kornr(i,1));
F1=x(kornr(i,1))-x(kornr(i,2));
G1=x(kornr(i,1))*(y(kornr(i,2))-y(kornr(i,1)))-y(kornr(i,1))*(x(kornr(i,2))-x(kornr(i,1)));
E2=y(kornr(i,4))-y(kornr(i,3));
F2=x(kornr(i,3))-x(kornr(i,4));
G2=x(kornr(i,3))*(y(kornr(i,4))-y(kornr(i,3)))-y(kornr(i,3))*(x(kornr(i,4))-x(kornr(i,3)));
xcord(i)=(-F2*G1+F1*G2)/(E2*F1-E1*F2);
ycord(i)=(-E2*G1+E1*G2)/(E1*F2-E2*F1);
DK(i,:)=[xcord(i),ycord(i)];
Juisteind(i)=i;
xcord1=xcord(i);
ycord1=ycord(i);
xcordp=xcord(i);
ycordp=ycord(i);
kount=kount+1;
i=i+1;
%
while not(i==num+1)
    if (kornr(i,1)==kornr(i-1,1))&&(kornr(i,2)==kornr(i-1,2))&&(kornr(i,3)==kornr(i-1,3))&&(kornr(i,4)==kornr(i-1,4))
        i=i+1;
    else
        if (kornr(i,1)==kornr(1,1))&&(kornr(i,2)==kornr(1,2))&&(kornr(i,3)==kornr(1,3))&&(kornr(i,4)==kornr(1,4))
            pp=find(not(Juisteind==0));
            ADK=DK(pp,:);
            empty=0;
            return
        else
            E1=y(kornr(i,2))-y(kornr(i,1));
            F1=x(kornr(i,1))-x(kornr(i,2));
            G1=x(kornr(i,1))*(y(kornr(i,2))-y(kornr(i,1)))-y(kornr(i,1))*(x(kornr(i,2))-x(kornr(i,1)));
            E2=y(kornr(i,4))-y(kornr(i,3));
            F2=x(kornr(i,3))-x(kornr(i,4));
            G2=x(kornr(i,3))*(y(kornr(i,4))-y(kornr(i,3)))-y(kornr(i,3))*(x(kornr(i,4))-x(kornr(i,3)));
            xcord(i)=(-F2*G1+F1*G2)/(E2*F1-E1*F2);
            ycord(i)=(-E2*G1+E1*G2)/(E1*F2-E2*F1);
            if ((abs(xcord(i)-xcordp)<eps)&&(abs(ycord(i)-ycordp)<eps))||((abs(xcord(i)-xcord1)<eps)&&(abs(ycord(i)-ycord1)<eps))
                i=i+1;
            else
                xcordp=xcord(i);
                ycordp=ycord(i);
                DK(i,:)=[xcord(i),ycord(i)];
                Juisteind(i)=i;
                kount=kount+1;
                i=i+1;
            end
        end
    end
end
%
% Delete all the empty spaces in the matrix DK. The result is ADK.
%
pp=find(not(Juisteind==0));
ADK=DK(pp,:);

%-------

function [tukmed]=halfmed(x,y,varargin)

% HALFMED is an algoritm that computes the Tukey median of a two-dimensional
% dataset. First, we have to check whether the data points are in general position.
% If not, a random number is added to each of the data points until the dataset
% comes in general position. All this is done in the m-file dithering. Then
% the deepest depth region is constructed. One can accelerate this search
% by giving the optional argument kstar, which is usually the maximal halfspace
% depth of the data points. This advantage is used in the functie bagp.m to
% quicken the computations. The Tukey median is the center of gravity of
% this deepest depth region.
%
% The HALFMED algoritm is described in:
%    Rousseeuw, P.J., Ruts, I. (1998),
%    "Constructing the bivariate Tukey median",
%    Statistica Sinica, 8, 827-839.
%
% Required input arguments:
%            x : vector containing the first coordinates of all the data
%                points
%            y : vector containing the second coordinates of all the data
%                points
%
% Optional input argument:
%        kstar : an integer between ceil(n/3) and floor(n/2)
%                One can also use the maximum halfspace depth of the data
%                points. (default = 0)
%
%
% I/O: result=halfmed(x,y,'kstar',0);
%  The name of the input arguments needs to be followed by their value.
%
% The output of halfmed is given by
%
%       result : the coordinates of the Tukey median
%
% This function is part of the Matlab Library for Robust Analysis,
% available at:
%              http://wis.kuleuven.be/stat/robust.html
%
% Last Update: 29/04/2005

xsum=0;
ysum=0;
tukmed=0;
n=length(x);
%
% Checking input
%
if length(x)==1
    error('x is not a vector')
elseif not(length(x)==length(y))
    error('The vectors x and y must have the same length.')
end
if sum(isnan(x))>=1 || sum(isnan(y))>=1
    error('Missing values are not allowed')
end
if sum(x==x(1))==n
    error('All data points ly on a vertical line.')
elseif sum(y==y(1))==n
    error('All data points ly on a horizontal line.')
else
    R=corrcoef(x,y);
    if abs(R(1,2))==1
        error('All data points are collineair.')
    end
end
%
% Looking for optional argument
%
if nargin>2
    if strcmp(varargin{1},'kstar')
        kstar=varargin{2};
        if length(kstar)>1
            error('kstar must be a number, not a vector.')
        end
    else
        error('Only kstar can be provided as an optional argument')
    end
else
    kstar=0;
end
if kstar>=floor(n/2)
    error('kstar must be smaller than floor(n/2) because the depth region corresponding with kstar is empty')
end
%
% Check whether the data are in general position. If not, add a very small random
% number to each of the data points.
%
[x,y, Index, angl, ind1,ind2]=dithering(x,y);
%
%Calculation of the Tukey median
%
if n<=3
    xsum=sum(x);
    ysum=sum(y);
    tukmed=[xsum/n,ysum/n];
    return
end
%
if kstar==0
    ib=ceil(n/3);
else
    ib=kstar;
end
ie=floor(n/2);
stay=1;
while stay==1
    le=ie-ib;
    if le<0
        le=0;
    end
    if le==0
        stay=0;
    end
    if stay==1
        [kou,dk,empty]=isodepth(x,y,ib+ceil(le/2));
        if empty==1
            ie=ib+ceil(le/2);
        end
        if empty==0
            ib=ib+ceil(le/2);
        end
        if le==1
            stay=0;
        end
    end
end
[kount,DK,empty]=isodepth(x,y,ib);
xsum=sum(DK(:,1));
if not(DK==0)
    ysum=sum(DK(:,2));
else
    ysum=0;
end
wx=DK(:,1);
if not(DK==0)
    wy=DK(:,2);
else
    wy=0;
end
%
% The maximal depth is now ib.
%
%
% Calculation of the center of gravity
%
som=0;
tukmed=0;
if kount>1
    wx=wx-(xsum/kount);
    wy=wy-(ysum/kount);
    for i=1:(kount-1)
        som=som+abs(wx(i)*wy(i+1)-wx(i+1)*wy(i));
        tukmed=tukmed+[(wx(i)+wx(i+1))*abs(wx(i)*wy(i+1)-wx(i+1)*wy(i)),(wy(i)+wy(i+1))*abs(wx(i)*wy(i+1)-wx(i+1)*wy(i))];
    end
    som=som+abs(wx(kount)*wy(1)-wx(1)*wy(kount));
    tukmed=tukmed+[(wx(kount)+wx(1))*abs(wx(kount)*wy(1)-wx(1)*wy(kount)),(wy(kount)+wy(1))*abs(wx(kount)*wy(1)-wx(1)*wy(kount))];
    tukmed=(tukmed/(3*som))+[(xsum/kount),(ysum/kount)];
else
    tukmed=[xsum,ysum];
end

%--------------

function result=bagp(x,y,sizesubset)

% BAGP is an algoritm that makes the necessary computations to construct
% the bagplot of a bivariate dataset. This function is especially used in
% bagplot.m. First, the data points are standardized. If the dataset is
% too large (n > 200), then a subset is taken. On the other hand, if the
% dataset is too small (n < 10), then no bag can be constructed. The Tukey
% median (see halfmed.m) is also calculated, so that we can center the data.
% Now, the bag lies between two consecutive depth regions. Therefore, one must
% search for the right depth k. Once the k-th and the (k-1)-th depth regions
% are constructed (see isodepth), the bag is the interpolation between the
% vertices of the 2 depth regions. Finally, all the data points are caracterized
% by their position compared to the bag with an integer 0,1,2 or 3.
%
% Required input arguments:
%            x : vector containing the first coordinates of all the data
%                points
%            y : vector containing the second coordinates of all the data
%                points
%
% I/O: result = bagp(x,y);
%
% The output of bagp is a structure containing
%
%        result.tukmed   : Tukey median
%        result.datatyp  : is 2 for observations in the bag, 1 for
%                          observations in the fence and 0 for outliers
%        result.interpol : The coordinates of the bag.
%        result.depth    : The halfspacedepth of each observation
%
% This function is part of the Matlab Library for Robust Analysis,
% available at:
%              http://wis.kuleuven.be/stat/robust.html
%
% Last Update: 29/04/2005

n=length(x);
eps=0.0000001;
nointer=0;
ntot=n;
%
% Checking input
%
if length(x)==1
    error('x is not a vector')
elseif not(length(x)==length(y))
    error('The vectors x and y must have the same length.')
end
if sum(isnan(x))>=1 || sum(isnan(y))>=1
    error('Missing values are not allowed')
end
if sum(x==x(1))==n
    error('All data points ly on a vertical line.')
elseif sum(y==y(1))==n
    error('All data points ly on a horizontal line.')
else
    R=corrcoef(x,y);
    if abs(R(1,2))==1
        error('All data points are collineair.')
    end
end
%
% Standardization of the data set
%
xmean=mean(x);
ymean=mean(y);
xdev=sqrt(var(x));
ydev=sqrt(var(y));
if xdev>eps
    x=(x-xmean)/xdev;
end
if ydev>eps
    y=(y-ymean)/ydev;
end
xoris=x;
yoris=y;
xori=x;
yori=y;
wx=x;
wy=y;
wx1=x;
wy1=y;
%
% If n is large, take a subset
%
nsub=sizesubset;
if n>nsub
    ntot=n;
    n=nsub;
    a=rdraw(n,ntot);
    x=wx1(a);
    y=wy1(a);
    wx=x;
    wy=y;
end
%
% Check whether the data is in general position. If not, add a very small random
% number to each of the data points.
%
[x,y, Index, angl, ind1,ind2]=dithering(x,y);
%
numdep(1:n)=0;
kstar=0;
for i=1:n
    hdep=halfspacedepth(x(i),y(i),x,y);
    if hdep>kstar
        kstar=hdep;
    end
    numdep(hdep)=numdep(hdep)+1;
end
%
% Calculation of the Tukey median
%
tukmed=halfmed(x,y,'kstar',kstar);
%
% Calculation of correct value of k
%
nc=floor(n/2);
j=kstar+1;
stay8=1;
while stay8==1
    stay8=0;
    j=j-1;
    if numdep(kstar)<=nc
        numdep(kstar)=numdep(kstar)+numdep(j-1);
        stay8=1;
    else
        k=j+1;
        numk1=numdep(kstar);
        numk=numk1-numdep(k-1);
        lambdanc=(nc-numk)/(numk1-numk);
    end
end
%
% Calculation of the vertices of Dk
%
[kount1,Dk1,empty1]=isodepth(x,y,k);
wx1=Dk1(:,1);
if Dk1==0
    wy1=0;
else
    wy1=Dk1(:,2);
end
%
% Calculation of the vertices of Dk-1
%
[kount2,Dk2,empty2]=isodepth(x,y,k-1);
wx2=Dk2(:,1);
if Dk2==0
    wy2=0;
else
    wy2=Dk2(:,2);
end
if corrcoef(wx2,wy2)>0.98
    error('Dk-1 is a line segment. Too many data points coincide.')
end
%
if n>=10
    wx1=wx1-tukmed(1);
    wy1=wy1-tukmed(2);
    wx2=wx2-tukmed(1);
    wy2=wy2-tukmed(2);
    x(1:n)=x(1:n)-tukmed(1);
    y(1:n)=y(1:n)-tukmed(2);
    xori(1:ntot)=xori(1:ntot)-tukmed(1);
    yori(1:ntot)=yori(1:ntot)-tukmed(2);
    %
    % Compute the angles of the data points.
    %
    numdattm=0;
    in=1:n;
    for i=1:n
        if (abs(x(i))<eps)&&(abs(y(i))<eps)
            numdattm=numdattm+1;
            dattm(numdattm)=i;
            angz(i)=1000;
        else
            dist=sqrt(x(i)^2+y(i)^2);
            xcord=x(i)/dist;
            ycord=y(i)/dist;
            if abs(xcord)>abs(ycord)
                if xcord>=0
                    angz(i)=asin(ycord);
                    if angz(i)<0
                        angz(i)=angz(i)+pi*2;
                    end
                else
                    angz(i)=pi-asin(ycord);
                end
            else
                if ycord>=0
                    angz(i)=acos(xcord);
                else
                    angz(i)=pi*2-acos(xcord);
                end
            end
            if angz(i)>=(2*pi-eps)
                angz(i)=0;
            end
        end
    end
    %
    % numdattm is the number of datapoints that are equal to the Tukey median
    %
    [angz,optel]=sort(angz);
    in=in(optel);
    n=n-numdattm;
    %
    % Compute the angles of the vertices of the depth region Dk.
    %
    if not(kount1==1)
        ind=1:kount1;
        for i=1:kount1
            dist=sqrt(wx1(i)^2+wy1(i)^2);
            if dist>eps
                xcord=wx1(i)/dist;
                ycord=wy1(i)/dist;
                if abs(xcord)>abs(ycord)
                    if xcord>=0
                        angy1(i)=asin(ycord);
                        if angy1(i)<0
                            angy1(i)=angy1(i)+pi*2;
                        end
                    else
                        angy1(i)=pi-asin(ycord);
                    end
                else
                    if ycord>=0
                        angy1(i)=acos(xcord);
                    else
                        angy1(i)=pi*2-acos(xcord);
                    end
                end
                if angy1(i)>=(pi*2-eps)
                    angy1(i)=0;
                end
            else
                nointer=1;
            end
        end
        [angy1,optel]=sort(angy1);
        ind=ind(optel);
    end
    %
    % Compute the angles of the vertices of the depth region Dk-1.
    %
    jnd=1:kount2;
    for i=1:kount2
        dist=sqrt(wx2(i)^2+wy2(i)^2);
        if dist>eps
            xcord=wx2(i)/dist;
            ycord=wy2(i)/dist;
            if abs(xcord)>abs(ycord)
                if xcord>=0
                    angy2(i)=asin(ycord);
                    if angy2(i)<0
                        angy2(i)=angy2(i)+pi*2;
                    end
                else
                    angy2(i)=pi-asin(ycord);
                end
            else
                if ycord>=0
                    angy2(i)=acos(xcord);
                else
                    angy2(i)=pi*2-acos(xcord);
                end
            end
            if angy2(i)>=(pi*2-eps)
                angy2(i)=0;
            end
        else
            nointer=1;
        end
    end
    [angy2,optel]=sort(angy2);
    jnd=jnd(optel);
    %
    % Calculation of arrays px and py for Dk
    %
    if not(kount1==1)
        jk=0;
        wx1(kount1+1)=wx1(1);
        wy1(kount1+1)=wy1(1);
        angy1(kount1+1)=angy1(1);
        ind(kount1+1)=ind(1);
        if angz(1)<angy1(1)
            j=kount1;
        end
        if angz(1)>=(angy1(1)-eps)
            j=1;
        end
        for i=1:n
            while (angz(i)>=(angy1(j+1)-eps))&&(jk==0)
                j=j+1;
                if j==kount1
                    jk=1;
                end
                if j==(kount1+1)
                    j=1;
                end
            end
            if (abs(wx1(ind(j+1))-wx1(ind(j)))>eps)&&(abs(x(in(i)))>eps)
                dist=y(in(i))/x(in(i))-(wy1(ind(j+1))-wy1(ind(j)))/(wx1(ind(j+1))-wx1(ind(j)));
                px(i,1)=(wy1(ind(j))*wx1(ind(j+1))-wx1(ind(j))*wy1(ind(j+1)))/(wx1(ind(j+1))-wx1(ind(j)));
                px(i,1)=px(i,1)/dist;
                py(i,1)=px(i,1)*y(in(i))/x(in(i));
            else
                if abs(wx1(ind(j+1))-wx1(ind(j)))<=eps
                    px(i,1)=wx1(ind(j));
                    py(i,1)=px(i,1)*y(in(i))/x(in(i));
                else
                    px(i,1)=0;
                    py(i,1)=((wy1(ind(j))-wy1(ind(j+1)))*wx1(ind(j))/(wx1(ind(j+1))-wx1(ind(j))))+wy1(ind(j));
                end
            end
        end
    end
    %
    % Calculation of arrays px and py for Dk-1
    %
    jk=0;
    wx2(kount2+1)=wx2(1);
    wy2(kount2+1)=wy2(1);
    angy2(kount2+1)=angy2(1);
    jnd(kount2+1)=jnd(1);
    if angz(1)<angy2(1)
        j=kount2;
    end
    if angz(1)>=(angy2(1)-eps)
        j=1;
    end
    for i=1:n
        while (angz(i)>=(angy2(j+1)-eps))&&(jk==0)
            j=j+1;
            if j==kount2
                jk=1;
            end
            if j==(kount2+1)
                j=1;
            end
        end
        if (abs(wx2(jnd(j+1))-wx2(jnd(j)))>eps)&&(abs(x(in(i)))>eps)
            dist=(y(in(i))/x(in(i)))-((wy2(jnd(j+1))-wy2(jnd(j)))/(wx2(jnd(j+1))-wx2(jnd(j))));
            px(i,2)=(wy2(jnd(j))*wx2(jnd(j+1))-wx2(jnd(j))*wy2(jnd(j+1)))/(wx2(jnd(j+1))-wx2(jnd(j)));
            px(i,2)=px(i,2)/dist;
            py(i,2)=px(i,2)*y(in(i))/x(in(i));
        else
            if abs(wx2(jnd(j+1))-wx2(jnd(j)))<=eps
                px(i,2)=wx2(jnd(j));
                py(i,2)=px(i,2)*y(in(i))/x(in(i));
            else
                px(i,2)=0;
                py(i,2)=((wy2(jnd(j))-wy2(jnd(j+1)))*wx2(jnd(j))/(wx2(jnd(j+1))-wx2(jnd(j))))+wy2(jnd(j));
            end
        end
    end
    if kount1==1
        px(:,1)=wx1(1);
        py(:,1)=wy1(1);
    end
    %
    % Mergesort of angy1 and angy2 to obtain gamma and calculation of arrays px
    % and py
    %
    if kount1==1
        gamma=angy2;
        px(:,3)=wx1;
        py(:,3)=wy1;
        px(:,4)=wx2;
        py(:,4)=wy2;
    else
        ja=1;
        jb=1;
        i=1;
        stay9=1;
        while stay9==1
            if angy1(ja)<=angy2(jb)
                if ja<=kount1
                    gamma(i)=angy1(ja);
                    px(i,3)=wx1(ind(ja));
                    py(i,3)=wy1(ind(ja));
                    if jb==1
                        wxjb1=wx2(jnd(kount2));
                        wyjb1=wy2(jnd(kount2));
                    else
                        wxjb1=wx2(jnd(jb-1));
                        wyjb1=wy2(jnd(jb-1));
                    end
                    if (abs(wx2(jnd(jb))-wxjb1)>eps)&&(abs(wx1(ind(ja)))>eps)
                        dist=(wy1(ind(ja))/wx1(ind(ja)))-((wy2(jnd(jb))-wyjb1)/(wx2(jnd(jb))-wxjb1));
                        px(i,4)=(wyjb1*wx2(jnd(jb))-wxjb1*wy2(jnd(jb)))/(wx2(jnd(jb))-wxjb1);
                        px(i,4)=px(i,4)/dist;
                        py(i,4)=px(i,4)*wy1(ind(ja))/wx1(ind(ja));
                    else
                        if abs(wx2(jnd(jb))-wxjb1)<=eps
                            px(i,4)=wxjb1;
                            py(i,4)=px(i,4)*wy1(ind(ja))/wx1(ind(ja));
                        else
                            px(i,4)=0;
                            py(i,4)=((wyjb1-wy2(jnd(jb)))*wxjb1/(wx2(jnd(jb))-wxjb1))+wyjb1;
                        end
                    end
                    ja=ja+1;
                    i=i+1;
                else
                    angy1(ja)=angy1(ja)+10;
                end
            else
                if jb<=kount2
                    gamma(i)=angy2(jb);
                    px(i,4)=wx2(jnd(jb));
                    py(i,4)=wy2(jnd(jb));
                    if ja==1
                        wxja1=wx1(ind(kount1));
                        wyja1=wy1(ind(kount1));
                    else
                        wxja1=wx1(ind(ja-1));
                        wyja1=wy1(ind(ja-1));
                    end
                    if (abs(wx1(ind(ja))-wxja1)>eps)&&(abs(wx2(jnd(jb)))>eps)
                        dist=(wy2(jnd(jb))/wx2(jnd(jb)))-((wy1(ind(ja))-wyja1)/(wx1(ind(ja))-wxja1));
                        px(i,3)=(wyja1*wx1(ind(ja))-wxja1*wy1(ind(ja)))/(wx1(ind(ja))-wxja1);
                        px(i,3)=px(i,3)/dist;
                        py(i,3)=px(i,3)*wy2(jnd(jb))/wx2(jnd(jb));
                    else
                        if abs(wx1(ind(ja))-wxja1)<=eps
                            px(i,3)=wxja1;
                            py(i,3)=px(i,3)*wy2(jnd(jb))/wx2(jnd(jb));
                        else
                            px(i,3)=0;
                            py(i,3)=((wyja1-wy1(ind(ja)))*wxja1/(wx1(ind(ja))-wxja1))+wyja1;
                        end
                    end
                    jb=jb+1;
                    i=i+1;
                else
                    angy2(jb)=angy2(jb)+10;
                end
            end
            if i<=(kount1+kount2)
                stay9=1;
            else
                stay9=0;
            end
        end
    end
    %
    % Interpolation of two arrays px and py
    %
    c=3;
    if nointer==1
        kount1=0;
    end
    if nointer==0
        wx2(1:(kount1+kount2))=lambdanc*px(1:(kount1+kount2),4)+(1-lambdanc)*px(1:(kount1+kount2),3);
        wy2(1:(kount1+kount2))=lambdanc*py(1:(kount1+kount2),4)+(1-lambdanc)*py(1:(kount1+kount2),3);
    end
    if xdev>eps
        xcord=(wx2(1:(kount1+kount2))+tukmed(1))*xdev+xmean;
    else
        xcord=(wx2(1:(kount1+kount2))+tukmed(1));
    end
    if ydev>eps
        ycord=(wy2(1:(kount1+kount2))+tukmed(2))*ydev+ymean;
    else
        ycord=wy2(1:(kount1+kount2))+tukmed(2);
    end
    interpol(1:(kount1+kount2),1)=xcord;
    interpol(1:(kount1+kount2),2)=ycord;
    %
    if nointer==1
        if xdev>eps
            xcord=(x(in(1:n))+tukmed(1))*xdev+xmean;
        else
            xcord=x(in(1:n))+tukmed(1);
        end
        if ydev>eps
            ycord=(y(in(1:n))+tukmed(2))*ydev+ymean;
        else
            ycord=y(in(1:n))+tukmed(2);
        end
        datatyp(1:n,1)=xcord;
        datatyp(1:n,2)=ycord;
        datatyp(1:n,3)=1;
        if xdev>eps
            xcord=(x(dattm(1:numdattm))+tukmed(1))*xdev+xmean;
        else
            xcord=x(dattm(1:numdattm))+tukmed(1);
        end
        if ydev>eps
            ycord=(y(dattm(1:numdattm))+tukmed(2))*ydev+ymean;
        else
            ycord=y(dattm(1:numdattm))+tukmed(2);
        end
        datatyp((n+1):(n+numdattm),1)=xcord;
        datatyp((n+1):(n+numdattm),2)=ycord;
        datatyp((n+1):(n+numdattm),3)=0;
        return
    end
    %
    % Repeat some calculations for the whole data set
    %
    if ntot>nsub
        n=ntot;
        x=xoris;
        y=yoris;
        %
        % Angles of the data points
        %
        x=x-tukmed(1);
        y=y-tukmed(2);
        numdattm=0;
        in=1:n;
        for i=1:n
            if (abs(x(i))<eps)&&(abs(y(i))<eps)
                angz(i)=1000;
                numdattm=numdattm+1;
                dattm(numdattm)=i;
            else
                dist=sqrt(x(i)^2+y(i)^2);
                xcord=x(i)/dist;
                ycord=y(i)/dist;
                if abs(xcord)>abs(ycord)
                    if xcord>=0
                        angz(i)=asin(ycord);
                        if angz(i)<0
                            angz(i)=angz(i)+pi*2;
                        end
                    else
                        angz(i)=pi-asin(ycord);
                    end
                else
                    if ycord>=0
                        angz(i)=acos(xcord);
                    else
                        angz(i)=pi*2-acos(xcord);
                    end
                end
                if angz(i)>=(pi*2-eps)
                    angz(i)=0;
                end
            end
        end
        [angz,optel]=sort(angz);
        in=in(optel);
        n=n-numdattm;
        %
        % Calculation of arrays px and py for B
        %
        jk=0;
        wx2(kount2+kount1+1)=wx2(1);
        wy2(kount2+kount1+1)=wy2(1);
        gamma(kount2+kount1+1)=gamma(1);
        if angz(1)<gamma(1)
            j=kount2+kount1;
        end
        if angz(1)>=(gamma(1)-eps)
            j=1;
        end
        for i=1:n
            while (angz(i)>=gamma(j+1)-eps)&&(jk==0)
                j=j+1;
                if j==(kount2+kount1)
                    jk=1;
                end
                if j==(kount2+kount1+1)
                    j=1;
                end
            end
            if (abs(wx2(j+1)-wx2(j))>eps)&&(abs(x(in(i)))>eps)
                dist=(y(in(i))/x(in(i)))-((wy2(j+1)-wy2(j))/(wx2(j+1)-wx2(j)));
                px(i,1)=(wy2(j)*wx2(j+1)-wx2(j)*wy2(j+1))/(wx2(j+1)-wx2(j));
                px(i,1)=px(i,1)/dist;
                py(i,1)=px(i,1)*y(in(i))/x(in(i));
            else
                if abs(wx2(j+1)-wx2(j))<=eps
                    px(i,1)=wx2(j);
                    py(i,1)=px(i,1)*y(in(i))/x(in(i));
                else
                    px(i,1)=0;
                    py(i,1)=((wy2(j)-wy2(j+1))*wx2(j)/(wx2(j+1)-wx2(j)))+wy2(j);
                end
            end
        end
    end
    %
    % Decide the type of each data point
    %
    c=3;
    num=0;
    num1=0;
    num2=0;
    num3=0;
    lopen=1;
    i=1;
    while lopen==1
        if ntot<=nsub
            px(i,1)=lambdanc*px(i,2)+(1-lambdanc)*px(i,1);
            py(i,1)=lambdanc*py(i,2)+(1-lambdanc)*py(i,1);
        end
        if (px(i,1)^2+py(i,1)^2)>eps
            lambda(i)=sqrt(x(in(i))^2+y(in(i))^2)/sqrt(px(i,1)^2+py(i,1)^2);
        else
            if (abs(x(in(i)))<eps) && (abs(y(in(i)))<eps)
                lambda(i)=0;
            else
                lambda(i)=c+1;
            end
        end
        if lambda(i)<eps
            num=num+1;
            typ(i)=0;
        elseif lambda(i)<=(1+eps)
            num1=num1+1;
            typ(i)=1;
        elseif lambda(i)<=(c+eps)
            num2=num2+1;
            typ(i)=2;
        elseif lambda(i)>c
            num3=num3+1;
            typ(i)=3;
        end
        i=i+1;
        if i==(n+1)
            lopen=0;
        end
    end
    if xdev>eps
        xcord=(x(in(1:n))+tukmed(1))*xdev+xmean;
    else
        xcord=x(in(1:n))+tukmed(1);
    end
    if ydev>eps
        ycord=(y(in(1:n))+tukmed(2))*ydev+ymean;
    else
        ycord=y(in(1:n))+tukmed(2);
    end
    datatyp(1:n,1)=xcord;
    datatyp(1:n,2)=ycord;
    typ=typ.';
    datatyp(1:n,3)=typ;
    %
    if not(numdattm==0)
        num=num+numdattm;
        typ(n+(1:numdattm))=0;
        lambda(n+(1:numdattm))=0;
        px(n+(1:numdattm),1)=x(dattm);
        py(n+(1:numdattm),1)=y(dattm);
        if xdev>eps
            xcord=(x(dattm)+tukmed(1))*xdev+xmean;
        else
            xcord=x(dattm)+tukmed(1);
        end
        if ydev>eps
            ycord=(y(dattm)+tukmed(2))*ydev+ymean;
        else
            ycord=y(dattm)+tukmed(2);
        end
        datatyp(n+(1:numdattm),1)=xcord;
        datatyp(n+(1:numdattm),2)=ycord;
        datatyp(n+(1:numdattm),3)=typ(n+(1:numdattm));
    end
else
    datatyp((1:n),1)=(x+tukmed(1))*xdev+xmean;
    datatyp((1:n),2)=(y+tukmed(2))*ydev+ymean;
    datatyp((1:n),3)=0;
    interpol=0;
end
%
if xdev>eps
    tukmed(1)=tukmed(1)*xdev+xmean;
end
if ydev>eps
    tukmed(2)=tukmed(2)*ydev+ymean;
end
result=struct('tukmed',tukmed,'datatyp',datatyp,'interpol',interpol,'depth',hdep);



%-----------------------------
function [x,y,Index,angl,ind1,ind2]=dithering(x,y)
%DITHERING is used to check whether the data points are in general position.
% If not, then a very small random number is added to each of
% the data points. Also the angles formed by pairs of data points are
% computed here. This file is used in the functions isodepth.m, halfmed.m
% and bagp.m.
%
% Required input arguments:
%            x : vector containing the first coordinates of all the data
%                points
%            y : vector containing the second coordinates of all the data
%                points

if not(length(x)==length(y))
    error('The vectors x and y must have the same length.')
end
if sum(isnan(x))>=1 || sum(isnan(y))>=1
    error('Missing values are not allowed')
end
n=length(x);
i=0;
blijf3=1;
dith=1;
xorig=x;
yorig=y;
while dith==1
    [xs,Index]=sort(x);
    ys=y(Index);
    dith=0;
    while blijf3==1
        blijf3=0;
        i=i+1;
        if (i+1)>n
            blijf3=2;
        else
            j=i+1;
            if not(xs(i)==xs(j))
                blijf3=1;
            else
                if ys(i)==ys(j)
                    dith=1;
                    blijf3=0;
                    % two datapoints coincide
                else
                    if ys(i)<ys(j)
                        even=Index(j);
                        Index(j)=Index(i);
                        Index(i)=even;
                    end
                    if (j+1)<=n
                        next=j+1;
                        if xs(i)==xs(next)
                            dith=1;
                            blijf3=0;
                            % three data points are collinear
                        else
                            blijf3=1;
                        end
                    end
                end
            end
        end
        if dith==1
            fac=1000000;
            ran=randn(n,2);
            x=xorig+ran(:,1)/fac;
            wx=x;
            y=yorig+ran(:,2)/fac;
            wy=y;
        else
            x=x;
            y=y;
        end
    end
    %
    if dith==0
        %
        % Compute all the angles formed by pairs of data points
        %
        m=0;
        for i=1:n
            rest=(i+1):n;
            spec=intersect(find(x(i)==x),rest);
            hoek(spec+m-i)=pi/2;
            spec=intersect(find(not(x(i)==x)),rest);
            hoek(spec+m-i)=atan((y(i)-y(spec))./(x(i)-x(spec)));
            p=m;
            m=m+n-i;
            if sum(hoek<=0)>=1
                spec=find(hoek<=0);
                hoek(spec)=hoek(spec)+pi;
            end
            ind1((p+1):m)=i;
            ind2((p+1):m)=rest;
        end
        %
        % Sort all the angles and permute ind1 and ind2 in the same way.
        %
        [angl,In]=sort(hoek);
        %
        ind1=ind1(In);
        ind2=ind2(In);
        %
        % Test wether any three datapoints are collinear.
        %
        ppp=diff(angl);
        k=find(ppp==0);
        %
        if sum(ind1(k)==ind1(k+1))>=1
            % There are 3 or more datapoints collineair.
            dith=1;
            fac=100000000;
            ran=randn(n,2);
            x=xorig+ran(:,1)/fac;
            wx=x;
            y=yorig+ran(:,2)/fac;
            wy=y;
        else
            x=x;
            y=y;
        end

    end
end


