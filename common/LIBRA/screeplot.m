function screeplot(vec,attrib,lev)

%SCREEPLOT draws the eigenvalues of the covariance matrix of the data in decreasing order.
% It can be used to decide how many latent variables to retain in a PCA analysis
% by looking at the kink in the curve. 
% 
% Reference: J.A. Jackson (1991), "A User's Guide to Principal Components",
%  Wiley Series in Probability and Mathematical Statistics, New-York.
%
% I/O: screeplot(vec,attrib,lev)
%
% Required input arguments:
%       vec : eigenvalues 
%    attrib : class identifier
% 
% Optional input argument: 
%    lev : 0 (default) or 1 to include a plot of the logarithm of the eigenvalues 
%
% Created 28 July 2000 by Sabine Verboven
% Last Update: 20/06/2003

if nargin < 3 
   lev=0;
end
if nargin < 2
    error('no class identifier is given.')
end
%%%%%%%MAIN%%%%%%%
set(gcf,'Name', 'Scree plot','NumberTitle', 'off');
cla
hold on
ymin=max([0,min(vec)]);
ymax=max(vec);
ymarg=0.06*(ymax-ymin);
ymin=ymin-ymarg;
ymax=ymax+ymarg;
plot(vec,'o-'); 
xlabel('Index')
ylabel('Eigenvalue')
if length(vec)==1
    xlim([0,2]); 
    ymarg=0.5;
    ymin=ymin-ymarg;
    ymax=ymin+(2*ymarg);
else 
    xlim([1,length(vec)]);
end
ylim([ymin,ymax]);
if lev==1
    set(gcf,'Name', 'LEV plot','NumberTitle', 'off');
    hold on
    plot(log(vec),'b.-'); 
    xlabel('Number of LV')
    ylabel('Log(eigenvalue)')
end

box on
title(attrib)
hold off