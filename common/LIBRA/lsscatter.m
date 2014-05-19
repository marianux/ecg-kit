function lsscatter(x,y,fitted,attrib)

%LSSCATTER makes a scatter plot with regression (LTS/LS) line
%
% Required input arguments:
%      x : predictor variabele (without missing values)
%      y : response variabele
% fitted : fitted values corresponding with the regression
% attrib : string identifying the used method =  'LS', 'LTS'
%
% I/O: lsscatter(x,y,fitted,attrib)
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Nele Smets on : 26/11/2003
% Last Update: 29/01/2008


set(gcf,'Name', 'Scatter plot', 'NumberTitle', 'off');
[n,p]=size(x);
if p~=1 
    disp(['Scatter plot with regression line ',...
            'is only available for bivariate data'])
else
    x1 = x(1);
    xn = x(n);
    y1 = fitted(1);
    yn = fitted(n);
    dx = xn - x1;
    dy = yn - y1;
    slope = dy./dx;
    centerx = (x1 + xn)/2;
    centery = (y1 + yn)/2;
    maxx = max(x);
    minx = min(x);
    maxy = centery + slope.*(maxx - centerx);
    miny = centery - slope.*(centerx - minx);
    mx = [minx; maxx];
    my = [miny; maxy];  
    plot(x,y,'o'); 
    hold on
    xrange=0.1*(abs(maxx-abs(minx)));
    yrange=0.2*(abs(maxy-abs(miny)+5));
    xlim([minx-xrange maxx+xrange])
    ylim([min(y)-yrange max(y)+yrange])
    plot(mx,my,'-');
    title(attrib)
    hold off
end

 