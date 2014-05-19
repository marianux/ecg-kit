function regresdiagplot3d(sdist,odist,rdist,cutoffsd,cutoffod,cutoffrd,k,class,multi,labsd,labod,labresd,labels)

%REGRESDIAGPLOT3D is a 3D-outlier map which visualizes the orthogonal distance, 
% the score distance and the residual distance calculated in a PCR or PLSR analysis.
%
% I/O: regresdiagplot3D(sdist,odist,rdist,cutoffsd,cutoffrd,k,class,labsd,labod,labresd,labels)
%
% Example: out=rpcr(X,Y);
%          regresdiagplot3d(out.sd,out.od,out.resd,out.cutoff.sd,cutoff.od,...
%          out.cutoff.resd,out.k,out.class,1,out.labsd,out.labod,out.labresd,0)
%
% Uses function: putlabel
%
% Created on by S.Verboven
% Last revision: 09/04/2004
% 

%INITIALIZATION%
if nargin<13
    labels=0;
end
if nargin<10
    labsd=3;
    labod=3;
    labresd=3;
    labels=0;
end
if nargin==10 
    labsd=3
    labels=0;
    labresd=3;
    labod=3;
end
if nargin==11
    labels=0;
    labresd=3;
    labod=3;
end
if nargin==12
    labels=0;
    labresd=3;
end
if nargin<9
    error('A required input variable is missing!')
end
    
%all LTS-analysis in RPCR are intercept included!!!
% if ask==2 %multivariate analysis
%     %residual distances
%    
% else  %univariate analysis
%     %standardized residuals
%     cutoffz=2.5; 
% end
cutoffxx=cutoffsd;
cutoffxy=cutoffsd; 
cutoffyy=cutoffod;
cutoffyx=cutoffod;
cutoffz=cutoffrd;
x=sdist;
y=odist;
z=rdist; 


%%%%%%%MAIN FUNCTION%%%%%%%
set(gcf,'Name', '3D-Outlier map (regression)', 'NumberTitle', 'off')%,'Renderer','OpenGL');

%%%%%%%%Odist=0 not yet included in this standalone!!!!!!!!! included in
%%%%%%%%makeplot function!!!!
plot3(x,y,z,'ko','markerfacecolor',[0.75 0.75 0.75])

hold on
axhandle=gca;
ylen=get(axhandle, 'Ylim');
xlen=get(axhandle,'Xlim');
zlen=get(axhandle,'Zlim');
xrange=xlen(2)-xlen(1);
upLimx=max(cutoffxx,xlen(2))+xrange*0.1;
lowLimx=xlen(1)-xrange*0.1;
yrange=ylen(2)-ylen(1);
upLimy=ylen(2)+yrange*0.1;
lowLimy=ylen(1)-yrange*0.1;
zrange=zlen(2)-zlen(1);
upLimz=max(cutoffz,zlen(2))+zrange*0.1;
if cutoffz==2.5
    lowLimz=min(-2.5,zlen(1))-zrange*0.1;
else
    lowLimz=min(cutoffz,zlen(1))-zrange*0.1;
end
    
%axis square;
set(gca, 'Xlim',[lowLimx upLimx],'Ylim',[lowLimy upLimy],'Zlim',[lowLimz,upLimz])
hold on

xlabel('Score distance')
ylabel('Orthogonal distance')
if cutoffz~=2.5
    zlabel('Residual distance')
else
    zlabel('Standardized Residual')
end
grid on


%in XY-space
%red plane "vertical" on x axis
oppy=[upLimy:-0.01:lowLimy];
n=length(oppy);
h=(upLimz-lowLimz)/n;
oppz=[lowLimz:h:upLimz];

X1=cutoffxx*ones(n);
Y1=repmat(oppy,n,1);
Z1=repmat([oppz(1:n-1)'; upLimz],1,n);
surf(X1,Y1,Z1,'edgecolor','none','facecolor','r')
alpha(.3)


%blue "horizontal" planes orthogonal on z-axis
%in ZX and ZY-space
oppx=[lowLimx:0.05:upLimx];
n=length(oppx);
X2=repmat(oppx,n,1);
h=(upLimy-lowLimy)/n;
oppy=[lowLimy:h:upLimy];
Y2=repmat(oppy(1:n)',1,n);
Z2=cutoffz*ones(n);
surf(X2,Y2,Z2,'edgecolor','none','facecolor','b')
alpha(0.3)

% in XY-space
%green plane "vertical" on y axis
oppx=[upLimx:-0.01:lowLimx];
n=length(oppx);
h=(upLimz-lowLimz)/n;
oppz=[lowLimz:h:upLimz];
X1=repmat(oppx,n,1);
Y1=cutoffyy*ones(n);
Z1=repmat([oppz(1:n-1)'; upLimz],1,n);
surf(X1,Y1,Z1,'edgecolor','none','facecolor','g')
alpha(0.3)

if multi==0 %univariate case
    surf(X2,Y2,-Z2,'edgecolor','none','facecolor','b')
    alpha(0.3)
end

if labels~=0
    putlabel(x,y,labels,z,labels)
else
    plotnumbers(x,y,labsd,labod,5,z,labresd)
end
hold off