function result=adjustedboxplot(x,varargin)
%
%ADJUSTEDBOXPLOT produces an adjusted box and whisker plot with one box for each column
% of X.  
% Typical for this boxplot are its skewness-adjusted whiskers, which are based on
% the medcouple, a robust measure of skewness. At skewed data, the original boxplot
% typically marks too many regular observations as outliers. The adjusted boxplot on the other
% hand makes a better distinction between regular observations and real outliers.
% 
% The ADJUSTED BOXPLOT is constructed as follows:
%   - a line is put at the height of the sample median.
%   - the box is drawn from the first to the third quartile.
%   - At right skewed data (MC >= 0), all points outside the interval 
%     [Q1 - 1.5*e^(a*MC)*IQR ; Q3 + 1.5*e^(b*MC)*IQR] are marked as outliers, 
%     where Q1 and Q3 denote the first and third quartile respectively, 
%     IQR stands for the interquartile range and MC is an abbreviation
%     for the medcouple (see mc.m).  
%     At left skewed data (MC < 0), the interval [Q1 - 1.5*e^(-b*MC)*IQR ;
%     Q3 + 1.5*e^(-a*MC)*IQR] is used, because of symmetry reasons.  
%   - finally, the whiskers are drawn, going from the ends of the box to the 
%     most remote points that are no outliers.  
% Note that the standard boxplot has the same box, but other whiskers. In that case, 
% all points outside the interval [Q1 - 1.5*IQR ; Q3 + 1.5*IQR] are marked as outliers. 
% The adjusted boxplot is thus similar to the original boxplot at symmetric distributions
% where MC=0.
%
% The skewness-adjusted boxplot is introduced in:
%     Hubert, M. and Vandervieren, E. (2007), 
%     "An adjusted boxplot for skewed distributions", 
%     Computational Statistics and Data Analysis, to appear.
%
% For the up-to-date reference, please consult the website:
%    http://wis.kuleuven.be/stat/robust.html
%
%
% Required input arguments:
%                 x : a data matrix; for each column the adjusted boxplot is
%                     drawn.  When also a grouping vector is given, x must be a
%                     vector.  
%
% Optional input arguments:
%                a : number, defining the whisker of the
%                    adjusted boxplot.  The default is to use a = -4.  
%                b : number, defining the whisker of the
%                    adjusted boxplot.  The default is to use b = 3.  
%                    This means that for right skewed data the interval 
%                    [Q1 - 1.5*e^(-4*MC)*IQR ; Q3 + 1.5*e^(3*MC)*IQR] and
%                    for left skewed data the interval
%                    [Q1 - 1.5*e^(-3*MC)*IQR ; Q3 + 1.5*e^(4*MC)*IQR] is
%                    used.  
%       groupvalid : Grouping variable defined as a vector, string matrix, 
%                    or cell array of strings.  Groupvalid can also be a 
%                    cell array of several grouping variables (such as 
%                    {G1 G2 G3}) to group the values in x by each unique 
%                    combination of grouping variable values.  
%          classic : If equal to 1, dotted lines are drawn at the height 
%                    of the original whiskers. If equal to 0, only the
%                    adjusted boxplot is plotted.  (default is 0).
%           symbol : Symbol and color to use for all outliers (default is 'r+').
%             vert : Box orientation, value 1 for a vertical boxplot (default) 
%                     or value 0 for a horizontal boxplot.
%           labels : Character array or cell array of strings containing
%                    labels for each column of x, or each group in G.
%           colors : A string or a three-column matrix of box colors.  Each
%                    box (outline, median line, and whiskers) is drawn in the
%                    corresponding color.  Default is to draw all boxes with
%                    blue outline, red median, and black whiskers.  Colors are
%                    recycled if mecessary.
%           widths : A numeric vector or scalar of box widths.  Default is
%                    0.5, or slightly smaller for fewer than three boxes.
%                    Widths are recycled if necessary.
%        positions : A numeric vector of box positions.  Default is 1:n.
%       grouporder : When G is given, a character array or cell array of
%                    group names, specifying the ordering of the groups in
%                    G.  Ignored when G is not given.
%
%
% I/O: result=adjustedboxplot(x,'a',-4,'b',3,'groupvalid',[],'classic',0,'symbol','r+',...
%             'orientation',1,'labels',[],'colors',[],'widths',1.5,'positions',[],'grouporder',[]);
%  The user should only give the input arguments that have to change their default value.
%  The name of the input arguments needs to be followed by their value.
%  The order of the input arguments is of no importance.
%
% ADJUSTEDBOXPLOT calls ADJUSTEDBOXUTIL to do the actual plotting.
%
% Examples: ADJUSTED BOXPLOT of car mileage grouped by country
%      load carsmall
%      out=adjustedboxplot(MPG,'groupvalid',Origin,'classic',1)
%      out=adjustedboxplot(MPG,'groupvalid',Origin,'symbol','b*','orientation',0)
%      out=adjustedboxplot(MPG,'groupvalid',Origin,'widths',0.75,'positions',[1 3 4 7 8 10], ...
%                   'grouporder',{'France' 'Germany' 'Italy' 'Japan' 'Sweden' 'USA'},'colors','kbrgym')
%
% The output is a structure containing the following fields:
%              result.a : numeric value, used in the definition of the outlier cutoffs 
%                         of the adjusted boxplot. 
%              result.b : numeric value, used in the definition of the outlier cutoffs 
%                         of the adjusted boxplot. 
%     result.groupvalid : If there is a grouping vector, this vector
%                         contains the assigned group numbers for the observations in vector x.
%        result.classic : value 1 or 0, indicating if dotted lines have been plot 
%                         at the heigth of the standard whiskers. 
%         result.symbol : symbol and color that has been used for all outliers.  
%    result.orientation : If equal to 1, the boxes are vertically oriented.  
%                         If equal to 0, a horizontal orientation has been used.  
%         result.labels : character array or cell array of strings, containing
%                         labels for each column of X, or each group in G.  
%         result.colors : a string or three-column matrix of box colors.
%         result.widths : a numeric vector of scalar of box widths.
%      result.positions : a numeric vector, containing the box positions.  
%     result.grouporder : If there is a grouping vector, this character array 
%                         or cell array of group names specifies the ordering of the groups.  
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at:
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Ellen Vandervieren
% Created on: 07/04/2005
% Last Update: 12/12/2006
%

if nargin<1
    error('Input argument ''x'' is undefined.')
end

whissw = 0; % don't plot whisker inside the box.

if isvector(x)
   % Might have one box, or might have a grouping variable. n will be properly
   % set later for the latter.
   x = x(:);
   n = 1; %
else
   % Have a data matrix, use as many boxes as columns.
   n = size(x,2);
end

% Assigning default-values
counter=1;
default=struct('a',-4,'b',3,'groupvalid',[],'classic',0,'symbol','r+','orientation',1,'labels',[],...
   'colors',[],'positions',[],'widths',0.5,'grouporder',[]);
% colors: default is blue box, red median, black whiskers
% positions: default is 1:n
% widths: default is 0.5, smaller for three or fewer boxes
% grouporder: default is 1:n
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

a=result.a;
b=result.b;
group=result.groupvalid;
classic=result.classic;
symbol=result.symbol;
orientation=result.orientation;
labels=result.labels;
colors=result.colors;
positions=result.positions;
widths=result.widths;
grouporder=result.grouporder;

% a and b must be numeric scalars.
if isempty(a)
    a = -4;
elseif ~isscalar(a) || ~isnumeric(a)
   error('LIBRA:adjustedboxplot:BadA',...
         'The ''a'' parameter value must be a numeric scalar.');
end

if isempty(b)
    b = 3;
elseif ~isscalar(b) || ~isnumeric(b)
   error('LIBRA:adjustedboxplot:BadB',...
         'The ''b'' parameter value must be a numeric scalar.');
end

% When group is non-empty, x must be a vector.  
if (~isempty(group) && ~isvector(x))
      error('LIBRA:adjustedboxplot:VectorRequired',...
            'x must be a vector when there is a grouping variable.');
end

% Classic must be equal to 0 or 1.
if isempty(classic)
   classic = 0;
elseif ~isscalar(classic) || ~ismember(classic,0:1)
   error('LIBRA:adjustedboxplot:InvalidClassic','Invalid value for ''classic'' parameter');
end

% Convert wordy inputs to internal codes
if isempty(orientation)
   orientation = 1;
elseif ischar(orientation)
   orientation = strmatch(orientation,{'horizontal' 'vertical'}) - 1;
end
if isempty(orientation) || ~isscalar(orientation) || ~ismember(orientation,0:1)   
   error('LIBRA:adjustedboxplot:InvalidOrientation',...
         'Invalid value for ''orientation'' parameter');
end

% Deal with grouping variable before processing more inputs
if ~isempty(group)
    if orientation, sep = '\n';
    else
        sep = ',';
    end
    [group,glabel,gname,multiline] = mgrp2idx(group,size(x,1),sep);
    n = size(gname,1);
    if numel(group) ~= numel(x)
        error('LIBRA:adjustedboxplot:InputSizeMismatch',...
            'X and G must have the same length.');
    end
else
    multiline = false;
end

% Reorder the groups if necessary
if ~isempty(group) && ~isempty(grouporder)
   if iscellstr(grouporder) || ischar(grouporder)
      % If we have a grouping vector, grouporder may be a list of group names.
      if ischar(grouporder), grouporder = cellstr(grouporder); end
      [dum,grouporder] = ismember(grouporder(:),glabel);
      % Must be a permutation of the group names
      if ~isequal(sort(grouporder),(1:n)')
         error('LIBRA:adjustedboxplot:BadGrouporder', ...
               'The ''grouporder'' parameter value must contain all the unique group names in G.');
      end
   else
      error('LIBRA:adjustedboxplot:BadGrouporder', ...
            'The ''grouporder'' parameter value must be a character array or a cell array of strings.');
   end
   group = order(group);
   glabel = glabel(grouporder);
   gname = gname(grouporder,:);
end

% Process the rest of the inputs

if isempty(labels)
   if ~isempty(group)
      labels = glabel;
   end
else
   if ~(iscellstr(labels) && numel(labels)==n) && ...
      ~(ischar(labels) && size(labels,1)==n)
      % Must have one label for each box
      error('LIBRA:adjustedboxplot:BadLabels','Incorrect number of box labels.');
   end
   if ischar(labels), labels = cellstr(labels); end
   multiline = false;
end
dfltLabs = (isempty(labels) && isempty(group)); % box labels are just column numbers

if isempty(widths)
   widths = repmat(min(0.15*n,0.5),n,1);
elseif ~isvector(widths) || ~isnumeric(widths) || any(widths<=0)
   error('LIBRA:adjustedboxplot:BadWidths', ...
         'The ''widths'' parameter value must be a numeric vector of positive values.');
elseif length(widths) < n
   % Recycle the widths if necessary.
   widths = repmat(widths(:),ceil(n/length(widths)),1);
end

if isempty(colors)
   % Empty colors tells adjustedboxutil to use defaults.
   colors = char(zeros(n,0));
elseif ischar(colors) && isvector(colors)
   colors = colors(:); % color spec string, make it a column
elseif isnumeric(colors) && (ndims(colors)==2) && (size(colors,2)==3)
   % RGB matrix, that's ok
else
   error('LIBRA:adjustedboxplot:BadColors',...
         'The ''colors'' parameter value must be a string or a three-column numeric matrix.');
end
if size(colors,1) < n
   % Recycle the colors if necessary.
   colors = repmat(colors,ceil(n/size(colors,1)),1);
end

if isempty(positions)
   positions = 1:n;
elseif ~isvector(positions) || ~isnumeric(positions)
   error('LIBRA:adjustedboxplot:BadPositions', ...
         'The ''positions'' parameter value must be a numeric vector.');
elseif length(positions) ~= n
   % Must have one position for each box
   error('LIBRA:adjustedboxplot:BadPositions', ...
         'The ''positions'' parameter value must have one element for each box.');
else
   if isempty(group) && isempty(labels)
      % If we have matrix data and the positions are not 1:n, we need to
      % force the default 1:n tick labels.
      labels = cellstr(num2str((1:n)'));
   end
end

%
% Done processing inputs
%

notch = 0;
whis = 1.5;

% Put at least the widest box or half narrowest spacing at each margin
if n > 1
    wmax = max(max(widths), 0.5*min(diff(positions)));
else
    wmax = 0.5;
end
xlims = [min(positions)-wmax, max(positions)+wmax];

ymin = nanmin(x(:));
ymax = nanmax(x(:));
if ymax > ymin
   dy = (ymax-ymin)/20;
else
   dy = 0.5;  % no data range, just use a y axis range of 1
end
ylims = [(ymin-dy) (ymax+dy)];

% Scale axis for vertical or horizontal boxes.
newplot
oldstate = get(gca,'NextPlot');
set(gca,'NextPlot','add','Box','on');
set(gcf,'Name', 'Adjusted boxplot', 'NumberTitle', 'off');
if orientation
    axis([xlims ylims]);
    set(gca,'XTick',positions);
    ylabel(gca,'Values');
    if dfltLabs, xlabel(gca, 'Column Number'); end
else
    axis([ylims xlims]);
    set(gca,'YTick',positions);
    xlabel(gca,'Values');
    if dfltLabs, ylabel(gca,'Column Number'); end
end
if nargout>0
   hout = [];
end

xvisible = NaN(size(x));
notnans = ~isnan(x);
for i= 1:n
   if ~isempty(group)
      thisgrp = find((group==i) & notnans);
   else
      thisgrp = find(notnans(:,i)) + (i-1)*size(x,1);
   end
   [outliers,hh] = adjustedboxutil(x(thisgrp),a,b,classic,notch,positions(i),widths(i), ...
                                 colors(i,:),symbol,orientation,whis,whissw);
   outliers = thisgrp(outliers);
   xvisible(outliers) = x(outliers);
   if nargout>0
      hout = [hout; hh(:)];
   end
end

if ~isempty(labels)
   if multiline && orientation
      % Turn off tick labels and axis label
      set(gca, 'XTickLabel','');
      setappdata(gca,'NLines',size(gname,2));
      xlabel(gca,'');
      ylim = get(gca, 'YLim');
      
      % Place multi-line text approximately where tick labels belong
      ypos = repmat(ylim(1),size(positions));
      text(positions,ypos,labels,'HorizontalAlignment','center', ...
                             'VerticalAlignment','top', 'UserData','xtick');

      % Resize function will position text more accurately
      set(gcf, 'ResizeFcn', @resizefcn, ...
               'Interruptible','off', 'PaperPositionMode','auto');
      resizefcn(gcf);
   elseif orientation
      set(gca, 'XTickLabel',labels);
   else
      set(gca, 'YTickLabel',labels);
   end
end
set(gca,'NextPlot',oldstate);

% Store information for gname function
set(gca, 'UserData', {'adjustedboxplot' xvisible group orientation});
hold off

%Output structure
result=struct('a',{a},'b',{b},'groupvalid',{group},'classic',{classic}, ...
    'symbol',{symbol},'orientation',{orientation},'labels',{labels},'colors',{colors},...
    'widths',{widths},'positions',{positions},'grouporder',{grouporder});

%=============================================================================

function [outlier,hout] = adjustedboxutil(x,a,b,classic,notch,lb,lf,clr,symbol,orientation,whis,whissw)
%ADJUSTEDBOXUTIL Produces a single adjusted boxplot.

% define the median and the quantiles
pctiles = prctile(x,[25;50;75]);
q1 = pctiles(1,:);
med = pctiles(2,:);
q3 = pctiles(3,:);

% find the extreme values (to determine where whiskers appear)
medc = mc(x);
if medc>=0
    vloadj = q1-whis*exp(a*medc)*(q3-q1);  %Lower cutoff value for the adjusted boxplot.  
    loadj = min(x(x>=vloadj));
    vhiadj = q3+whis*exp(b*medc)*(q3-q1);  %Upper cutoff value for the adjusted boxplot.  
    upadj = max(x(x<=vhiadj));
else
    vloadj = q1-whis*exp(-b*medc)*(q3-q1);  %Lower cutoff value for the adjusted boxplot.  
    loadj = min(x(x>=vloadj));
    vhiadj = q3+whis*exp(-a*medc)*(q3-q1);  %Upper cutoff value for the adjusted boxplot.  
    upadj = max(x(x<=vhiadj));
end

if (isempty(loadj)), loadj = q1; end
if (isempty(upadj)), upadj = q3; end

if (isequal(classic,1)),
    vloorig = q1-whis*(q3-q1);     %Lower cutoff value for the original boxplot. 
    loorig = min(x(x>=vloorig));
    if (isempty(loorig)), loorig = q1; end
end

if (isequal(classic,1)),
    vhiorig = q3+whis*(q3-q1);     %Upper cutoff value for the original boxplot.      
    uporig = max(x(x<=vhiorig));
    if (isempty(uporig)), uporig = q3; end
end

x1 = repmat(lb,1,2);
x2 = x1+[-0.25*lf,0.25*lf];
outlier = x<loadj | x > upadj;
yy = x(outlier);

xx = repmat(lb,1,length(yy));
lbp = lb + 0.5*lf;
lbm = lb - 0.5*lf;

if whissw == 0
   upadj = max(upadj,q3);
   loadj = min(loadj,q1);
   if (isequal(classic,1)),
       uporig = max(uporig,q3);  
       loorig = min(loorig,q1);
   end
end

% Set up (X,Y) data for notches if desired.
if ~notch
    xx2 = [lbm lbp lbp lbm lbm];
    yy2 = [q3 q3 q1 q1 q3];
    xx3 = [lbm lbp];
else
    n1 = med + 1.57*(q3-q1)/sqrt(length(x));
    n2 = med - 1.57*(q3-q1)/sqrt(length(x));
    if n1>q3, n1 = q3; end
    if n2<q1, n2 = q1; end
    lnm = lb-0.25*lf;
    lnp = lb+0.25*lf;
    xx2 = [lnm lbm lbm lbp lbp lnp lbp lbp lbm lbm lnm];
    yy2 = [med n1 q3 q3 n1 med n2 q1 q1 n2 med];
    xx3 = [lnm lnp];
end
yy3 = [med med];

    
% Determine if the boxes are vertical or horizontal.
% The difference is the choice of x and y in the plot command.
if orientation
    if (isequal(classic,1)),
            hout = plot(x1,[q3 upadj],'k--', x1,[loadj q1],'k--',...
                   x2,[loadj loadj],'k-',...
                   x2,[upadj upadj],'k-', xx2,yy2,'b-',xx3,yy3,'r-',xx,yy,symbol,...
                   x2,repmat(loorig,1,length(x2)),'k:',x2,repmat(uporig,1,length(x2)),'k:');
            title('Adjusted and original boxplot');  
    else
            hout = plot(x1,[q3 upadj],'k--', x1,[loadj q1],'k--',...
                x2,[upadj upadj],'k-', x2,[loadj loadj],'k-', ...
                xx2,yy2,'b-', xx3,yy3,'r-', xx,yy,symbol);
            title('Adjusted boxplot');
    end
else
    if (isequal(classic,1)),
            hout = plot([q3 upadj],x1,'k--', [loadj q1],x1,'k--',...
                   [loadj loadj],x2,'k-', ...
                   [upadj upadj],x2,'k-',yy2,xx2,'b-',yy3,xx3,'r-',yy,xx,symbol,... 
                   x2,repmat(loorig,1,length(x2)),'k:',x2,repmat(uporig,1,length(x2)),'k:');
            title('Adjusted and original boxplot');   
    else
            hout = plot([q3 upadj],x1,'k--', [loadj q1],x1,'k--',...
                [upadj upadj],x2,'k-', [loadj loadj],x2,'k-', ...
                yy2,xx2,'b-', yy3,xx3,'r-', yy,xx,symbol);
            title('Adjusted boxplot');
    end    
end

% If there's a color given, show everything in that color.  If the outlier
% symbol has a color in it, leave those alone.
if ~isempty(clr)
   if any(ismember(symbol,'bgrcmyk'))
      set(hout(1:6),'Color',clr);
   else
      set(hout,'Color',clr);
   end
end

set(hout(1),'Tag','Upper Whisker');
set(hout(2),'Tag','Lower Whisker');
set(hout(3),'Tag','Upper Adjacent Value');
set(hout(4),'Tag','Lower Adjacent Value');
set(hout(5),'Tag','Box');
set(hout(6),'Tag','Median');
if length(hout)>=7
   set(hout(7),'Tag','Outliers');
else
   hout(7) = NaN;
end


%=============================================================================

function resizefcn(f)%,dum)

% Adjust figure layout to make sure labels remain visible
h = findobj(f, 'UserData','xtick');
if (isempty(h))
   set(f, 'ResizeFcn', '');
   return;
end
ax = get(f, 'CurrentAxes');
nlines = getappdata(ax, 'NLines');

% Position the axes so that the fake X tick labels have room to display
set(ax, 'Units', 'characters');
p = get(ax, 'Position');
ptop = p(2) + p(4);
if (p(4) < nlines+1.5)
   p(2) = ptop/2;
else
   p(2) = nlines + 1;
end
p(4) = ptop - p(2);
set(ax, 'Position', p);
set(ax, 'Units', 'normalized');

% Position the labels at the proper place
xl = get(ax, 'XLabel');
set(xl, 'Units', 'data');
p = get(xl, 'Position');
ylim = get(ax, 'YLim');
p2 = (p(2)+ylim(1))/2;
for j=1:length(h)
   p = get(h(j), 'Position') ;
   p(2) = p2;
   set(h(j), 'Position', p);
end

%==========================================================================
function [ogroup,glabel,gname,multigroup] = mgrp2idx(group,rows,sep)
%MGRP2IDX Convert multiple grouping variables to index vector
%   [OGROUP,GLABEL,GNAME,MULTIGROUP] = MGRP2IDX(GROUP,ROWS) takes
%   the inputs GROUP, ROWS, and SEP.  GROUP is a grouping variable (numeric
%   vector, string matrix, or cell array of strings) or a cell array
%   of grouping variables.  ROWS is the number of observations.
%   SEP is a separator for the grouping variable values.
%
%   The output OGROUP is a vector of group indices.  GLABEL is a cell
%   array of group labels, each label consisting of the values of the
%   various grouping variables separated by the characters in SEP.
%   GNAME is a cell array containing one column per grouping variable
%   and one row for each distinct combination of grouping variable
%   values.  MULTIGROUP is 1 if there are multiple grouping variables
%   or 0 if there are not.

%   Tom Lane, 12-17-99
%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 1.4.2.1 $  $Date: 2004/01/24 09:36:20 $

multigroup = (iscell(group) & size(group,1)==1);
if (~multigroup)
   [ogroup,gname] = grp2idx(group);
   glabel = gname;
else
   % Group according to each distinct combination of grouping variables
   ngrps = size(group,2);
   grpmat = zeros(rows,ngrps);
   namemat = cell(1,ngrps);
   
   % Get integer codes and names for each grouping variable
   for j=1:ngrps
      [g,gn] = grp2idx(group{1,j});
      grpmat(:,j) = g;
      namemat{1,j} = gn;
   end
   
   % Find all unique combinations
   [urows,ui,uj] = unique(grpmat,'rows');
   
   % Create a cell array, one col for each grouping variable value
   % and one row for each observation
   ogroup = uj;
   gname = cell(size(urows));
   for j=1:ngrps
      gn = namemat{1,j};
      gname(:,j) = gn(urows(:,j));
   end
   
   % Create another cell array of multi-line texts to use as labels
   glabel = cell(size(gname,1),1);
   if (nargin > 2)
      nl = sprintf(sep);
   else
      nl = sprintf('\n');
   end
   fmt = sprintf('%%s%s',nl);
   lnl = length(fmt)-3;        % one less than the length of nl
   for j=1:length(glabel)
      gn = sprintf(fmt, gname{j,:});
      gn(end-lnl:end) = [];
      glabel{j,1} = gn;
   end
end
