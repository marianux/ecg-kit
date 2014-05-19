function simcaplot(result); 

%SIMCAPLOT plots a scatter plot with the boundaries defined by the SIMCA method. 
% It is based on the results from a simca analysis (see rsimca.m or csimca.m).
% 
% For technical reasons, 6 different groups can be plotted (with different symbols).
% In case there are more groups, please adapt lines 15--17.
%
% I/O: simcaplot(result)
%
% Written by K. Vanden Branden on 01/08/04
% Last update: 21/09/2004

p = size(result.x,2);
markings={'bx';'ro';'kdiamond';'y+';'g.';'m*'};
linecolor={'b-','r-','k-','y-','g-','m-'};
colorellipse={'b','r','k','y','g','m'};
if p>3
    error('The dimension of the dataset is larger than 3.')
elseif length(result.avemisclasprob)>6
    error(['Only 6 groups can be drawn with different symbols. Please adapt the code',... 
            ' of simcaplot.m if you want to draw more groups.'])
else
    legstr=[];
    nClass = length(result.groupmisclasprob);
    for iClass = 1:nClass
        out = result.pca{iClass};
        groupi = result.x(find(result.group == iClass),:);
        %legstr=[legstr; sprintf(['Group',num2str(iClass)])]; 
        if out.k == 1
            [out.lbord,out.rbord]=confreg(out);
            bord=[out.lbord,out.rbord];
        elseif out.k == 2
            bord=ellipse([0 0],diag(out.L(1:2)));
        elseif out.k == 3
            bord=ellips3D(out.M, out.P*diag(out.L)*out.P');
        end
        if p == 2
            xlabel('x');ylabel('y');
            if size(bord) == [1,4]
                plot(groupi(:,1),groupi(:,2),markings{iClass});hold on;
                plot([bord(1,1),bord(1,3)],[bord(1,2),bord(1,4)],linecolor{iClass});
            elseif size(bord)==[802,2]
                plot(groupi(:,1),groupi(:,2),markings{iClass});hold on;
                newbord=bord*out.P'+repmat(out.M,802,1);
                line(newbord(:,1),newbord(:,2));
                
            end
        elseif p == 3
            xlabel('x');ylabel('y');zlabel('z')
            if size(bord) == [1,6]
                plot3(groupi(:,1),groupi(:,2),groupi(:,3),markings{iClass});hold on;grid on;
                plot3([bord(1,1),bord(1,4)],[bord(1,2),bord(1,5)],[bord(1,3),bord(1,6)],linecolor{iClass})
            elseif size(bord)==[802,2]
                plot3(groupi(:,1),groupi(:,2),groupi(:,3),markings{iClass});hold on;grid on;
                newbord=bord*out.P'+repmat(out.M,802,1);
                plot3(newbord(:,1),newbord(:,2),newbord(:,3),linecolor{iClass});
            else 
                plot3(groupi(:,1),groupi(:,2),groupi(:,3),markings{iClass});hold on;grid on;
                mesh(bord.x,bord.y1,bord.z1,'Edgecolor',colorellipse{iClass});alpha(0.2);
            end
        end
    end
    title(result.class)
    %legend(legstr,0)
    hold off
end


%------------------
function [lbord,rbord,outl]=confreg(result);

n1=size(result.T,1);p=size(result.P,1);
rbord=result.M+sqrt(result.L)*result.cutoff.sd*result.P';
lbord=result.M-sqrt(result.L)*result.cutoff.sd*result.P';

%----------------------
function coord=ellipse(mean,covar)

% Determines the coordinates of some points that lie on the 97.5 % tolerance ellipse.

deter=covar(1,1)*covar(2,2)-covar(1,2)^2;
ylimit=sqrt(7.37776*covar(2,2));
y=-ylimit:0.005*ylimit:ylimit;
sqtdi=sqrt(deter*(ylimit^2-y.^2))/covar(2,2);
sqtdi([1,end])=0;
b=mean(1)+covar(1,2)/covar(2,2)*y;
x1=b-sqtdi;
x2=b+sqtdi;
y=mean(2)+y;
coord=[x1,x2([end:-1:1]);y,y([end:-1:1])]';

%----------------
function out = ellips3D(mu, sigma)

%This function determines the 3D ellipsoid.

dist = sqrt(chi2inv(0.975,3));

A = inv(sigma);
[ev,ew] = eig(A);
A = ev'*A*ev;

c1 = 0.1;
zlimit = dist/sqrt(A(3,3));
z = -zlimit:(c1*zlimit):zlimit;
c2 = (1/c1*2 + 1)*2; % is length(z)*2
c3 = c2/2; % is  length(z)
z1 = ones(1,c2)*z(1);
x1 = zeros(c3,c3);
x2 = zeros(c3,c3);
y = zeros(c3,c3);
y1 = zeros(c3,c2);

for i = 2:(c3-1)
    z1 = [z1 ones(1,c2)*z(i)];
    D = dist^2 - A(3,3)*z(i)^2;
    ylimit = sqrt(D)/sqrt(A(2,2));
    y(i,:) = -ylimit:(c1*ylimit):ylimit;
    y1(i,:) = [y(i,:) ylimit:(-c1*ylimit):-ylimit];
    
    for j = 2:(c3-1)
        x1(i,j) = sqrt((D - A(2,2)*y(i,j)^2)/A(1,1));
        x2(i,j) = -sqrt((D - A(2,2)*y(i,j)^2)/A(1,1));
    end
end

x = [x1,x2];
z1 = [z1 ones(1,c2)*z(c3)];
z1 = vec2mat(z1,c2);
grotex = vec2mat([mat2vec(x) mat2vec(y1) mat2vec(z1)],c2*c3);
groottex = ev*grotex;

groottex(1,:) = groottex(1,:) + mu(1);
groottex(2,:) = groottex(2,:) + mu(2);
groottex(3,:) = groottex(3,:) + mu(3);

x = vec2mat(groottex(1,:),c2);
y1 = vec2mat(groottex(2,:),c2);
z1 = vec2mat(groottex(3,:),c2);
out.x=x;
out.y1=y1;
out.z1=z1;


%-----------------
function vec = mat2vec(mat)

nkolom = size(mat,2);
nrij = size(mat,1);

vec = [];

for i=1:nkolom
   hulpvec = wkeep(mat,[nrij,1],[1,1]);
   vec = [vec hulpvec'];
   mat = mat(:,2:(nkolom - (i-1)));
end

%----------------------

function mat = vec2mat(vec,ncol)

nrow = length(vec)/ncol;

for i = 1:nrow
   mat(i,:)  = wkeep(vec,ncol,'l');
   vec = vec((ncol + 1):length(vec));
end