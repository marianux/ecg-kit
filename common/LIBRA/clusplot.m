function [] = clusplot(x,ncluv)

%CLUSPLOT creates a bivariate plot visualizing a partition (clustering)
% of the data. All observations are represented by points in the plot,
% using principal components or multidimensional scaling. Around each
% cluster an ellipse is drawn.
%
%The algorithm is fully described in:
%   Pison, G., Struyf, A., and Rousseeuw, P.J. (1999),
%   "Displaying a clustering with CLUSPLOT",
%   Computational Statistics and Data Analysis,
%   30, 381-392.
%
% Required input arguments:
%   x : Data matrix (rows = observations, columns = variables)
%       All variables need to be numeric.
%       Dissimilarity matrix (if input is vector)
%   ncluv : Cluster membership for each observation (mostly
%           as been found by PAM, FANNY or CLARA)
%
% I/O:
%   result=clusplot(x,ncluv)
%
% Example (subtracted from the referenced paper)
%   load ruspini.mat
%   result=pam(ruspini,4,[4 4]);
%   clusplot(ruspini,result.ncluv);
%   clusplot(result.dys,result.ncluv);
%
% The output of CLUSPLOT is a figure containing the data
% (possibly plotted in a lower dimension) and ellipses
% surrounding the clusters.
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at:
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Guy Brys (May 2006)

close all
span=1;
xgrp=double(ncluv);

if (isempty(x))
    error('The input matrix must be non-empty');
elseif ((size(x,1)==1) | (size(x,2)==1))
    if (size(x,1)==1)
        x = x';
    end
    %Dissimilarities given
    if (sum(isnan(x))~=0)
        error('Missing values not allowed in the dissimilarity vector');
    end
    %We assume x to be a vector of dissimilarities
    n = (1+sqrt(1+8*size(x,1)))/2;
    x1 = zeros(n,n);
    indexcol = 0;
    for i=1:(n-1)
        indexcol = [indexcol ; (max(indexcol)-i+1)+repmat(i,(n-i),1)+[i:(n-1)]'];
    end
    indexcol=indexcol(2:length(indexcol))';
    x1(indexcol) = x;
    x1 = x1+x1';
    addc = sclprc(x1,length(ncluv),1);
    [evec,eiga] = eig(addc);
    eiga = eiga(1:2,1:2);
    evec = evec(:,1:2);
    points = -evec*sqrt(eiga);
    if (sum(eiga<0)>0)
        error('There are non-positive eigenvalues');
    end
    plotgrp(points(:,1),points(:,2),xgrp);
    set(gcf,'Position',[340 150 720 460]);
    title(sprintf('CLUSPLOT (the %d components explain %5.2f percent of the point variability)',2,max(0,min(100,100*sum(sum(eiga))/sum(diag(addc))))))
    hold on
    xplot = points
elseif (size(x,2)==2)
    plotgrp(x(:,1),x(:,2),xgrp);
    set(gcf,'Position',[340 150 720 460]);
    hold on
    xplot = x;
else
    %Real data given
    if ((sum((sum(isnan(x))==size(x,1)))>0) | (sum((sum(isnan(x'))==size(x,2)))>0))
        error('A variable or observation contains only missing values');
    end
    indexnan = find(isnan(x));
    for i=indexnan
        x(i) = median(x(isnan(x(:,ceil(i/size(x,1))))~=1,ceil(i/size(x,1))));
    end
    R = cpca(x,'k',2,'plots',0);
    R.T(:,2) = -R.T(:,2);
    plotgrp(R.T(:,1),R.T(:,2),xgrp);
    set(gcf,'Position',[340 150 720 460]);
    title(sprintf('CLUSPLOT (the %d components explain %5.2f percent of the point variability)',2,100*sum(R.L)/sum(var(x))))
    hold on
    xplot = R.T;
end

k = double(max(ncluv));

for i=1:k
    x1 = xplot(ncluv==i,:);
    if (rank(cov(x1))==2)
        if (span==1)
            x2 = [repmat(1,size(x1,1),1) x1];
            [clsqdi,clprob,clstop] = spannc(x2);
            [B(i,:) A(:,:,i)] = weightmecov_new(x1,max(0,clprob));
            D(i) = sqrt(weightmecov_new(clsqdi',max(0,clprob)));
        elseif (span==0)
            B(i,:) = mean(x1);
            A(:,:,i) = cov(x1);
            D(i) = sqrt(max(mahalanobis(x1,B(i,:),'cov',A(:,:,i))));
            D(i) = D(i)+0.01*D(i);
        end
    elseif (rank(cov(x1))~=2)
        if (span==1)
            D(i) = 1;
            if ((sum(x1(:,1)~=x1(1,1))~=0) | (sum(x1(:,2)~=x1(1,2))~=0))
                B(i,:) = [min(x1(:,1))+(max(x1(:,1))-min(x1(:,1)))/2 min(x1(:,2))+(max(x1(:,2))-min(x1(:,2)))/2];
                aa = sqrt((B(i,1)-min(x1(:,1)))^2+(B(i,2)-min(x1(:,2)))^2);
                if (sum(x1(:,1)~=x1(1,1))~=0)
                    [maxx1,ind1] = max(x1(:,1));
                    [minx1,ind2] = min(x1(:,1));
                    qq = atan((x1(ind1,2)-x1(ind2,2))/(maxx1-minx1));
                    bb = 0;
                else
                    qq = pi/2;
                    bb = 0;
                end
                A(:,:,i) = [cos(qq) sin(qq) ; -sin(qq) cos(qq)]*[aa^2 0 ; 0 bb^2]*[cos(qq) -sin(qq) ; sin(qq) cos(qq)];
            else
                aa = (max(x1(:,1))-min(x1(:,1)))/90;
                bb = (max(x1(:,2))-min(x1(:,2)))/70;
                aa = aa+(aa==0);
                bb = bb+(bb==0);
                A(:,:,i) = [aa^2 0 ; 0 bb^2];
                B(i,:) = x1(1,:);
            end
        else
            D(i) = 1;
            if (((max(x1(:,1))-min(x1(:,1)))>(max(x(:,1))-min(x(:,1)))/70) | ((max(x1(:,2))-min(x1(:,2)))>(max(x(:,2))-min(x(:,2)))/50))
                B(i,:) = [min(x1(:,1))+(max(x1(:,1))-min(x1(:,1)))/2 min(x1(:,2))+(max(x1(:,2))-min(x1(:,2)))/2];
                aa = sqrt((B(i,1)-min(x1(:,1)))^2+(B(i,2)-min(x1(:,2)))^2);
                aa = aa+0.05*aa;
                if ((max(x1(:,1))-min(x1(:,1)))>(max(x(:,1))-min(x(:,1)))/70)
                    [maxx1,ind1] = max(x1(:,1));
                    [minx1,ind2] = min(x1(:,1));
                    qq = atan((x1(ind1,2)-x1(ind2,2))/(maxx1-minx1));
                    if (min(x(:,2))==max(x(:,2)))
                        bb = 1;
                    else
                        if ((max(x1(:,2))-min(x1(:,2)))>(max(x(:,2))-min(x(:,2)))/50)
                            bb = (max(x1(:,2))-min(x1(:,2)))/10;
                        else
                            bb = (max(x(:,2))-min(x(:,2)))/40;
                        end
                    end
                else
                    if (min(x(:,1))==max(x(:,1)))
                        bb = 1;
                    else
                        bb = (max(x(:,1))-min(x(:,1)))/40;
                    end
                    qq = pi/2;
                end
                A(:,:,i) = [cos(qq) sin(qq) ; -sin(qq) cos(qq)]*[aa^2 0 ; 0 bb^2]*[cos(qq) -sin(qq) ; sin(qq) cos(qq)];
            else
                aa = (max(x(:,1))-min(x(:,1)))/90;
                bb = (max(x(:,2))-min(x(:,2)))/70;
                aa = aa+(aa==0);
                bb = bb+(bb==0);
                A(:,:,i) = [aa^2 0 ; 0 bb^2];
                B(i,:) = x1(1,:);
            end
        end
    end
    posedges = ellipse(A(:,:,i),B(i,:),D(i));
    patch(posedges(:,1),posedges(:,2),[0 0 0],'FaceColor','none')
end
hold off

%------------
%SUBFUNCTIONS

function [posedges] = ellipse(a,b,d)

deta = a(1,1)*a(2,2)-a(1,2)^2;
ylimit = sqrt(a(2,2))*d;
y = -ylimit:(0.01*ylimit):ylimit;
discr = sqrt((deta/(a(2,2)^2))*(a(2,2)*d^2-y.^2));
discr([1 length(discr)]) = 0;
z = b(1)+(a(1,2)/a(2,2))*y;
x1 = z-discr;
x2 = z+discr;
y = b(2)+y;
posedges = [x1' y' ; x2(length(x2):-1:1)' y(length(y):-1:1)'];

%---
function [] = plotgrp(x,y,ncluv)

symb = ['+' 'o' 's' 'd' 'x' '*' '.' 'v'];
k = max(ncluv);
if (k<=8)
    for i=1:k
        plot(x(ncluv==i),y(ncluv==i),symb(i));
        hold on
    end
else
    error('Too many groups');
end

%---
function [wmean,wcov]=weightmecov_new(data,weights)

[n,nvar] = size(data);
if (~(isempty(find(weights<0))))
    error('The weights are negative');
end
if size(weights,1)==1
    weights=weights';
end
wmean=sum(diag(weights/sum(weights))*data);

wcov=((data - repmat(wmean,n,1))'*diag(weights)*(data - repmat(wmean,n,1)))/(1-sum(weights.^2));
