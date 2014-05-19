function daplot(x,group,center,covar,classic,method)

%DAPLOT plots 97.5% tolerances ellipses of the bivariate data set x, which
% consists of several groups defined by the input argument 'group'. 
% Center and covar are estimates of the center and covariance matrix of each group, 
% obtained with a classical ('CDA') or a robust ('RDA') discriminant analysis. See cda.m and rda.m
% Method indicates whether linear ('linear') or quadratic ('quadratic')
% classification is used. 
%
% For technical reasons, 6 different groups can be plotted (with different symbols).
% In case there are more groups, please adapt lines 22 and 23.
%
% I/O: daplot(x,group,center,cov,'CDA','linear')
%
% Written by S. Verboven on 28/07/03
% Last update: 17/02/2004

set(gcf,'Name', 'Discriminant analysis', 'NumberTitle', 'off');
n = size(x,1);
[m,p] = size(center);
lev=unique(group);
marking={'bx';'ro';'kdiamond';'y+';'g.';'m*'};
linecolor={'b-','r-','k-','y-','g-','m-'};
if p~=2
    disp('Warning: Tolerance ellipses are only drawn for two-dimensional data sets.')
    return
elseif length(lev)>6
    error(['Only 6 groups can be drawn with different symbols. Please adapt the code',... 
    ' of daplot.m if you want to draw more groups.'])
else
    legstr=[];
    for i=1:m
        data{i}=x(group==lev(i),:);
        plot(data{i}(:,1),data{i}(:,2),marking{i});   
        hold on
        legstr=[legstr; sprintf(['Group',num2str(i)])];
    end 
    for i=1:m
        if strcmp(method,'linear')
            Z{i} = ellipse(center(i,:),covar); 
        elseif strcmp(method,'quadratic')
            Z{i} = ellipse(center(i,:),covar{i}); 
        end
        plot(Z{i}(:, 1), Z{i}(:, 2), linecolor{i})
        hold on
    end
end
if strmatch('CDA',classic,'exact') & strmatch('linear',method,'exact')
    title('Classical linear discriminant analysis')
elseif strmatch('CDA',classic,'exact') & strmatch('quadratic',method,'exact')
    title('Classical quadratic discriminant analysis')  
elseif strmatch('RDA',classic,'exact') & strmatch('linear',method,'exact')
    title('Robust linear discriminant analysis')
elseif strmatch('RDA',classic,'exact') & strmatch('quadratic',method,'exact')
    title('Robust quadratic discriminant analysis')
end
hold off
legend(legstr,0)

%--------------
function out = ellipse(loc,covar)

A=covar;
detA = A(1, 1)*A(2, 2) - A(1, 2)^2;
dist = sqrt(chi2inv(0.975, 2));
ylimit = sqrt(A(2, 2)) * dist;
y = - ylimit:0.01*ylimit:ylimit;
sqroot.discr = sqrt(detA/A(2, 2).^2 * (A(2, 2) * dist.^2 - y.^2));
sqroot.discr([1, length(sqroot.discr)]) = 0;
b = loc(1) + A(1, 2)/A(2, 2) * y;
x1 = b - sqroot.discr;
x2 = b + sqroot.discr;
y = loc(2) + y;
out= [[x1' y']; flipud([x2' y'])]; 

