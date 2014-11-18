function varargout=plot2svg(param1,id,pixelfiletype)
%  Matlab to SVG converter
%  Prelinary version supporting 3D plots as well
%
%  Usage: plot2svg(filename,graphic handle,pixelfiletype)
%                  optional     optional     optional
%         or
%
%         plot2svg(figuresize,graphic handle,pixelfiletype)
%                   optional     optional      optional
%
%         pixelfiletype = 'png' (default), 'jpg'
%
%  Juerg Schwizer 23-Oct-2005
%  See http://www.juergschwizer.de to get more informations
%
%  07.06.2005 - Bugfix axxindex (Index exceeds matrix dimensions)
%  19.09.2005 - Added possibility to select output format of pixel graphics
%  23.10.2005 - Bugfix cell array strings (added by Bill)
%               Handling of 'hggroups' and improved grouping of objects
%               Improved handling of pixel images (indexed and true color pictures)
%  23.10.2005 - Switched default pixelfromat to 'png'
%  07.11.2005 - Added handling of hidden axes for annotations (added by Bill)
%  03.12.2005 - Bugfix of viewBox to make Firefox 1.5 working
%  04.12.2005 - Improved handling of exponent values for log-plots
%               Improved markers
%  09.12.2005 - Bugfix '<' '>' '?' '"'
%  22.12.2005 - Implementation of preliminary 3D version
%               Clipping
%               Minor tick marks
%  22.01.2005 - Removed unused 'end'
%  29.10.2006 - Bugfix '°','±','µ','²','³','¼''½','¾','©''®'
%  17-04-2007 - Bugfix 'projection' in hggroup and hgtransform
%  27-01-2008 - Added Octave functionality (thanks to Jakob Malm)
%               Bugfixe cdatamapping (thanks to Tom)
%               Bugfix image data writing (thanks to Tom)
%               Patches includes now markers as well (needed for 'scatter'
%               plots (thanks to Phil)
%  04-02-2008 - Bugfix markers for Octave (thanks to Jakob Malm)
%  30-12-2008 - Bugfix image scaling and orientation
%               Bugfix correct backslash (thanks to Jason Merril)
%  20-06-2009 - Improvment of image handling (still some remaining issues)
%               Fix for -. line style (thanks to Ritesh Sood)
%  28-06-2009 - Improved depth sorting for patches and surface
%             - Bugfix patches
%             - Bugfix 3D axis handling

global PLOT2SVG_globals
global colorname
progversion='28-Jun-2009';
PLOT2SVG_globals.runningIdNumber = 0;
PLOT2SVG_globals.octave = false;
if nargout==1
    varargout={0};
end
disp(['   Matlab/Octave to SVG converter version ' progversion ', Juerg Schwizer (converter@juergschwizer.de).'])
matversion=version;
if exist('OCTAVE_VERSION','builtin')
    PLOT2SVG_globals.octave = true;
    disp('   Info: PLOT2SVG runs in Octave mode.')
else
    if str2num(matversion(1))<6 % Check for matlab version and print warning if matlab version lower than version 6.0 (R.12)
        disp('   Warning: Future versions may no more support older versions than MATLAB R12.')
    end
end
if nargout > 1
    error('Function returns only one return value.')
end
if nargin<2 % Check if handle was included into function call, otherwise take current figure
    id=gcf;
end
if nargin==0
    if PLOT2SVG_globals.octave
        error('PLOT2SVG in Octave mode does not yet support a file menu. File name is needed during function call.')
    else
        [filename, pathname] = uiputfile( {'*.svg', 'SVG File (*.svg)'},'Save Figure as SVG File');
        if ~( isequal( filename, 0) || isequal( pathname, 0))    
            % yes. add backslash to path (if not already there)
            pathname = addBackSlash( pathname); 
            % check, if extension is allrigth
            if ( ~strcmpi( getFileExtension( filename), '.svg'))
                filename = [ filename, '.svg'];
            end
            finalname=[pathname filename];
        else
            disp('   Cancel button was pressed.')
            return
        end
    end
else
    if isnumeric(param1)
        if PLOT2SVG_globals.octave
            error('PLOT2SVG in Octave mode does not yet support a file menu. File name is needed during function call.')
        else
            [filename, pathname] = uiputfile( {'*.svg', 'SVG File (*.svg)'},'Save Figure as SVG File');  
            if ~( isequal( filename, 0) || isequal( pathname, 0))    
                % yes. add backslash to path (if not already there)
                pathname = addBackSlash( pathname); 
                % check, if ectension is allrigth
                if ( ~strcmpi( getFileExtension( filename), '.svg'))
                    filename = [ filename, '.svg'];
                end
                finalname=[pathname filename];
            else
                disp('   Cancel button was pressed.')
                return
            end
        end
    else
        finalname=param1;   
    end
end
% needed to see annotation axes
originalShowHiddenHandles = get(0, 'ShowHiddenHandles');
set(0, 'ShowHiddenHandles', 'on');
originalFigureUnits=get(id,'Units');
set(id,'Units','pixels');   % All data in the svg-file is saved in pixels
paperpos=get(id,'Position');
if ( nargin > 0)
    if isnumeric(param1)
        paperpos(3)=param1(1);
        paperpos(4)=param1(2);
    end
end
if (nargin < 3)
    PLOT2SVG_globals.pixelfiletype = 'png';
else
    PLOT2SVG_globals.pixelfiletype = pixelfiletype;
end
cmap=get(id,'Colormap');
colorname='';
for i=1:size(cmap,1)
    colorname(i,:)=sprintf('%02x%02x%02x',fix(cmap(i,1)*255),fix(cmap(i,2)*255),fix(cmap(i,3)*255));
end

% Open SVG-file
[pathstr,name,ext] = fileparts(finalname);
%PLOT2SVG_globals.basefilename = fullfile(pathstr,name);
PLOT2SVG_globals.basefilepath = pathstr;
PLOT2SVG_globals.basefilename = name;
PLOT2SVG_globals.figurenumber = 1;
fid=fopen(finalname,'wt');   % Create a new text file
fprintf(fid,'<?xml version="1.0" standalone="no"?><!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n');    % Insert file header
fprintf(fid,'<svg preserveAspectRatio="xMinYMin meet" width="100%%" height="100%%" viewBox="0 0 %0.3f %0.3f" ',paperpos(3),paperpos(4));
fprintf(fid,'  version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">\n');
fprintf(fid,'  <desc>Matlab Figure Converted by PLOT2SVG written by Juerg Schwizer</desc>\n');
group=1;
groups=[];
axfound=0;
% Frame of figure
figcolor = searchcolor(id,get(id, 'Color'));
if (~ strcmp(figcolor, 'none'))
    % Draw rectangle in the background of the graphic frame to cover all
    % other graphic elements
    try % Octave does not have support for InvertHardcopy yet -- Jakob Malm
        if strcmp(get(id,'InvertHardcopy'),'on')
            fprintf(fid,'  <rect x="0" y="0" width="%0.3f" height="%0.3f" fill="#ffffff" stroke="none" />\n',paperpos(3),paperpos(4));
        else
            fprintf(fid,'  <rect x="0" y="0" width="%0.3f" height="%0.3f" fill="%s" stroke="none" />\n',paperpos(3),paperpos(4),figcolor);
        end
    catch
        fprintf(fid,'  <rect x="0" y="0" width="%0.3f" height="%0.3f" fill="%s" stroke="none" />\n',paperpos(3),paperpos(4),figcolor);
    end
end
% Search all axes
ax=get(id,'Children');
for j=length(ax):-1:1
    currenttype = get(ax(j),'Type');
    if strcmp(currenttype,'axes')
        group=group+1;
        groups=[groups group];
        group=axes2svg(fid,id,ax(j),group,paperpos);
        axfound=1;
    elseif strcmp(currenttype,'uicontrol')
        if strcmp(get(ax(j),'Visible'),'on')
            control2svg(fid,id,ax(j),group,paperpos);
            axfound=1;
        end
    elseif strcmp(currenttype, 'uicontextmenu') || ...
            strcmp(currenttype, 'uimenu') || ...
            strcmp(currenttype, 'hgjavacomponent') || ...
            strcmp(currenttype, 'uitoolbar')
        % ignore these types
    else
        disp(['   Warning: Unhandled main figure child type: ' currenttype]);
    end
end
fprintf(fid,'</svg>\n');
fclose(fid);    % close text file
if nargout==1
    varargout={0};
end
set(id,'Units',originalFigureUnits);
set(0, 'ShowHiddenHandles', originalShowHiddenHandles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS %%%%%
% Create axis frame and insert all children of this axis frame
function group=axes2svg(fid,id,ax,group,paperpos)
global colorname
global PLOT2SVG_globals
originalAxesUnits=get(ax,'Units');
set(ax,'Units','normalized');
axpos=get(ax,'Position');
faces =    [1 2 4 3; 2 4 8 6; 3 4 8 7; 1 2 6 5; 1 5 7 3; 5 6 8 7];
%           x-y    ; y-z    ; x-z    ; y-z    ; x-z    ; x-y
corners(:,:,1) = [1 1 2 3 4; 2 1 3 2 4];
corners(:,:,2) = [2 2 4 6 8; 3 2 6 4 8];
corners(:,:,3) = [1 3 4 7 8; 3 3 7 4 8];
corners(:,:,4) = [1 1 2 5 6; 3 1 5 2 6];
corners(:,:,5) = [2 1 3 5 7; 3 1 5 3 7];
corners(:,:,6) = [1 5 6 7 8; 2 5 7 6 8];
edge_neighbours = [2 3 5; 1 4 6; 4 1 7; 3 2 8; 6 7 1; 5 8 2; 8 5 3; 7 6 4];
edge_opposite = [8 7 6 5 4 3 2 1];
nomx = [0 1 0 1 0 1 0 1];
nomy = [0 0 1 1 0 0 1 1];
nomz = [0 0 0 0 1 1 1 1];
[projection,edges] = get_projection(ax,id);
x = (edges(1,:)*axpos(3)+axpos(1))*paperpos(3);
y = (1-(edges(2,:)*axpos(4)+axpos(2)))*paperpos(4);    
% Depth Sort of view box edges 
[edge_z,edge_index]=sort(edges(3,:));
most_back_edge_index = edge_index(1);
% Back faces are plot box faces that are behind the plot (as seen by the
% view point)
back_faces = find(any(faces == most_back_edge_index,2));
front_faces = find(all(faces ~= most_back_edge_index,2));
groupax=group;
axlimx=get(ax,'XLim');
axlimy=get(ax,'YLim');
axlimz=get(ax,'ZLim');
axlimxori=axlimx;
axlimyori=axlimy;
axlimzori=axlimz;
if strcmp(get(ax,'XScale'),'log')
    axlimx=log10(axlimx);
    axlimx(find(isinf(axlimx)))=0;
end
if strcmp(get(ax,'YScale'),'log')
    axlimy=log10(axlimy);
    axlimy(find(isinf(axlimy)))=0;
end
if strcmp(get(ax,'ZScale'),'log')
    axlimz=log10(axlimz);
    axlimz(find(isinf(axlimz)))=0;
end
if strcmp(get(ax,'XDir'),'reverse')
    axlimx = fliplr(axlimx);
end
if strcmp(get(ax,'YDir'),'reverse')
    axlimy = fliplr(axlimy);
end
if strcmp(get(ax,'ZDir'),'reverse')
    axlimz = fliplr(axlimz);
end
axlimori = [axlimxori(1) axlimyori(1) axlimzori(1) axlimxori(2)-axlimxori(1) axlimyori(2)-axlimyori(1) axlimzori(2)-axlimzori(1)];
fprintf(fid,'  <g>\n');
axIdString = createId;
fprintf(fid,'  <clipPath id="%s">\n',axIdString);
fprintf(fid,'    <rect x="%0.3f" y="%0.3f" width="%0.3f" height="%0.3f"/>\n',...
    min(x), min(y), max(x)-min(x), max(y)-min(y));
fprintf(fid,'  </clipPath>\n');
if strcmp(get(ax,'Visible'),'on')
    group=group+1;
    grouplabel=group;
    axxtick=get(ax,'XTick');
    axytick=get(ax,'YTick');
    axztick=get(ax,'ZTick');
    axlabelx=get(ax,'XTickLabel');
    axlabely=get(ax,'YTickLabel');
    axlabelz=get(ax,'ZTickLabel');
    if PLOT2SVG_globals.octave
        if projection.xyplane
            axlabelz = [];
        end
    end
    gridlinestyle=get(ax,'GridLineStyle');
    minor_gridlinestyle=get(ax,'MinorGridLineStyle');
    try % Octave does not have 'TickLength' yet. --Jakob Malm
        both_ticklength = get(ax,'TickLength');
    catch
        both_ticklength = [ 0.01 0.025 ];
    end
    if projection.xyplane
        ticklength = both_ticklength(1);
        xy_ratio = axpos(3)*paperpos(3)/ (axpos(4)*paperpos(4));
        if xy_ratio < 1
            tick_ratio = [1 1/xy_ratio 1];
        else
            tick_ratio = [xy_ratio 1 1];
        end
        if strcmp(get(ax,'TickDir'),'out')
            label_distance = -(0.02 + ticklength);
        else
            label_distance = -0.02;
        end
    else
        ticklength = both_ticklength(2);
        label_distance = -2*abs(ticklength);
        tick_ratio = [1 1 1];
    end
    linewidth = get(ax,'LineWidth');
    axxindex=find((axxtick >= axlimori(1)) & (axxtick <= (axlimori(1)+axlimori(4))));
    axyindex=find((axytick >= axlimori(2)) & (axytick <= (axlimori(2)+axlimori(5))));
    axzindex=find((axztick >= axlimori(3)) & (axztick <= (axlimori(3)+axlimori(6))));
    % remove sticks outside of the axes (-1 of legends)
    axxtick=axxtick(axxindex); 
    axytick=axytick(axyindex);
    axztick=axztick(axzindex);
    if length(axxtick) > 1
        minor_lin_sticks = (0.2:0.2:0.8)*(axxtick(2)-axxtick(1));
        minor_axxtick = [];
        for stick = [2*axxtick(1)-axxtick(2) axxtick]
            minor_axxtick = [minor_axxtick minor_lin_sticks + stick]; 
        end
        minor_axxtick = minor_axxtick(find(minor_axxtick > min(axlimx) & minor_axxtick < max(axlimx)));
    else
        minor_axxtick = [];
    end
    if length(axytick) > 1
        minor_lin_sticks = (0.2:0.2:0.8)*(axytick(2)-axytick(1));
        minor_axytick = [];
        for stick = [2*axytick(1)-axytick(2) axytick]
            minor_axytick = [minor_axytick minor_lin_sticks + stick]; 
        end
        minor_axytick = minor_axytick(find(minor_axytick > min(axlimy) & minor_axytick < max(axlimy)));
    else
        minor_axytick = [];
    end
    if length(axztick) > 1
        minor_lin_sticks = (0.2:0.2:0.8)*(axztick(2)-axztick(1));
        minor_axztick = [];
        for stick = [2*axztick(1)-axztick(2) axztick]
            minor_axztick = [minor_axztick minor_lin_sticks + stick]; 
        end
        minor_axztick = minor_axztick(find(minor_axztick > min(axlimz) & minor_axztick < max(axlimz)));
    else
        minor_axztick = [];
    end
    axxindex_inner = find((axxtick > axlimori(1)) & (axxtick < (axlimori(1)+axlimori(4))));
    axyindex_inner = find((axytick > axlimori(2)) & (axytick < (axlimori(2)+axlimori(5))));
    axzindex_inner = find((axztick > axlimori(3)) & (axztick < (axlimori(3)+axlimori(6))));
    minor_log_sticks = log10(0.2:0.1:0.9);
    if strcmp(get(ax,'TickDir'),'out')
        ticklength=-ticklength;
        valid_xsticks = 1:length(axxindex);
        valid_ysticks = 1:length(axyindex);
        valid_zsticks = 1:length(axzindex);
    else
        valid_xsticks = axxindex_inner;
        valid_ysticks = axyindex_inner;
        valid_zsticks = axzindex_inner;
    end
    if strcmp(get(ax,'XScale'),'log')
        axxtick = log10(get(ax,'XTick'));
        minor_axxtick = [];
        for stick = axxtick
            minor_axxtick = [minor_axxtick minor_log_sticks + stick]; 
        end
        minor_axxtick = minor_axxtick(find(minor_axxtick > min(axlimx) & minor_axxtick < max(axlimx)));
    end
    if strcmp(get(ax,'YScale'),'log')
        axytick=log10(get(ax,'YTick'));
        minor_axytick = [];
        for stick = axytick
            minor_axytick = [minor_axytick minor_log_sticks + stick]; 
        end
        minor_axytick = minor_axytick(find(minor_axytick > min(axlimy) & minor_axytick < max(axlimy)));
    end
    if strcmp(get(ax,'ZScale'),'log')
        axztick=log10(get(ax,'ZTick'));
        minor_axztick = [];
        for stick = axztick
            minor_axztick = [minor_axztick minor_log_sticks + stick]; 
        end
        minor_axztick = minor_axztick(find(minor_axztick > min(axlimz) & minor_axztick < max(axlimz)));
    end
    % Draw back faces 
    linewidth=get(ax,'LineWidth');
    if ~strcmp(get(ax,'Color'),'none')
        background_color = searchcolor(id,get(ax,'Color'));
        background_opacity = 1;
    else
        background_color = '#000000';
        background_opacity = 0;
    end
    for p=1:size(back_faces)
        patch2svg(fid, group, axpos, x(faces(back_faces(p),:)), y(faces(back_faces(p),:)), background_color, '-', linewidth, 'none', background_opacity, 1.0)
    end
    for pindex = 1:size(back_faces)
        p = back_faces(pindex);
        for k = 1:size(corners,1)
            switch corners(k,1,p)
                case 1 % x
                    % Draw x-grid
                    scolorname=searchcolor(id,get(ax,'XColor'));
                    if strcmp(get(ax,'XGrid'),'on')
                        if axlimx(1)~=axlimx(2)
                            xg_line_start = interp1([axlimx(1) axlimx(2)],[x(corners(k,2,p)) x(corners(k,3,p))],axxtick);
                            yg_line_start = interp1([axlimx(1) axlimx(2)],[y(corners(k,2,p)) y(corners(k,3,p))],axxtick);
                            xg_line_end = interp1([axlimx(1) axlimx(2)],[x(corners(k,4,p)) x(corners(k,5,p))],axxtick);
                            yg_line_end = interp1([axlimx(1) axlimx(2)],[y(corners(k,4,p)) y(corners(k,5,p))],axxtick);
                            for i = axxindex_inner
                                line2svg(fid,grouplabel,axpos,[xg_line_start(i) xg_line_end(i)],[yg_line_start(i) yg_line_end(i)],scolorname,gridlinestyle,linewidth)
                            end
                            if strcmp(get(ax,'XTickMode'),'auto') && strcmp(get(ax,'XMinorGrid'),'on') && ~isempty(minor_axxtick)
                                xg_line_start = interp1([axlimx(1) axlimx(2)],[x(corners(k,2,p)) x(corners(k,3,p))],minor_axxtick);
                                yg_line_start = interp1([axlimx(1) axlimx(2)],[y(corners(k,2,p)) y(corners(k,3,p))],minor_axxtick);
                                xg_line_end = interp1([axlimx(1) axlimx(2)],[x(corners(k,4,p)) x(corners(k,5,p))],minor_axxtick);
                                yg_line_end = interp1([axlimx(1) axlimx(2)],[y(corners(k,4,p)) y(corners(k,5,p))],minor_axxtick);
                                for i = 1:length(xg_line_start)
                                    line2svg(fid,grouplabel,axpos,[xg_line_start(i) xg_line_end(i)],[yg_line_start(i) yg_line_end(i)],scolorname,minor_gridlinestyle,linewidth)
                                end  
                            end
                        end
                    end
                    if projection.xyplane == false
                        if strcmp(get(ax,'Box'),'on')
                            line2svg(fid,grouplabel,axpos,[x(corners(k,2,p)) x(corners(k,3,p))],[y(corners(k,2,p)) y(corners(k,3,p))],scolorname,'-',linewidth);
                            line2svg(fid,grouplabel,axpos,[x(corners(k,4,p)) x(corners(k,5,p))],[y(corners(k,4,p)) y(corners(k,5,p))],scolorname,'-',linewidth);
                        else
                            if strcmp(get(ax,'XGrid'),'on')
                                line2svg(fid,grouplabel,axpos,[x(corners(k,2,p)) x(corners(k,3,p))],[y(corners(k,2,p)) y(corners(k,3,p))],scolorname,gridlinestyle,linewidth);
                                line2svg(fid,grouplabel,axpos,[x(corners(k,4,p)) x(corners(k,5,p))],[y(corners(k,4,p)) y(corners(k,5,p))],scolorname,gridlinestyle,linewidth);
                            end
                        end
                    end
                case 2 % y
                    % Draw y-grid
                    scolorname=searchcolor(id,get(ax,'YColor'));
                    if strcmp(get(ax,'YGrid'),'on')
                        if axlimy(1)~=axlimy(2)
                            xg_line_start = interp1([axlimy(1) axlimy(2)],[x(corners(k,2,p)) x(corners(k,3,p))],axytick);
                            yg_line_start = interp1([axlimy(1) axlimy(2)],[y(corners(k,2,p)) y(corners(k,3,p))],axytick);
                            xg_line_end = interp1([axlimy(1) axlimy(2)],[x(corners(k,4,p)) x(corners(k,5,p))],axytick);
                            yg_line_end = interp1([axlimy(1) axlimy(2)],[y(corners(k,4,p)) y(corners(k,5,p))],axytick);
                            for i = axyindex_inner
                                line2svg(fid,grouplabel,axpos,[xg_line_start(i) xg_line_end(i)],[yg_line_start(i) yg_line_end(i)],scolorname,gridlinestyle,linewidth)
                            end
                            if strcmp(get(ax,'YTickMode'),'auto') && strcmp(get(ax,'YMinorGrid'),'on') && ~isempty(minor_axytick)
                                xg_line_start = interp1([axlimy(1) axlimy(2)],[x(corners(k,2,p)) x(corners(k,3,p))],minor_axytick);
                                yg_line_start = interp1([axlimy(1) axlimy(2)],[y(corners(k,2,p)) y(corners(k,3,p))],minor_axytick);
                                xg_line_end = interp1([axlimy(1) axlimy(2)],[x(corners(k,4,p)) x(corners(k,5,p))],minor_axytick);
                                yg_line_end = interp1([axlimy(1) axlimy(2)],[y(corners(k,4,p)) y(corners(k,5,p))],minor_axytick);
                                for i = 1:length(xg_line_start)
                                    line2svg(fid,grouplabel,axpos,[xg_line_start(i) xg_line_end(i)],[yg_line_start(i) yg_line_end(i)],scolorname,minor_gridlinestyle,linewidth)
                                end  
                            end
                        end
                    end
                    if projection.xyplane == false
                        if strcmp(get(ax,'Box'),'on')
                            line2svg(fid,grouplabel,axpos,[x(corners(k,2,p)) x(corners(k,3,p))],[y(corners(k,2,p)) y(corners(k,3,p))],scolorname,'-',linewidth);
                            line2svg(fid,grouplabel,axpos,[x(corners(k,4,p)) x(corners(k,5,p))],[y(corners(k,4,p)) y(corners(k,5,p))],scolorname,'-',linewidth);
                        else
                            if strcmp(get(ax,'YGrid'),'on')
                                line2svg(fid,grouplabel,axpos,[x(corners(k,2,p)) x(corners(k,3,p))],[y(corners(k,2,p)) y(corners(k,3,p))],scolorname,gridlinestyle,linewidth);
                                line2svg(fid,grouplabel,axpos,[x(corners(k,4,p)) x(corners(k,5,p))],[y(corners(k,4,p)) y(corners(k,5,p))],scolorname,gridlinestyle,linewidth);
                            end
                        end
                    end
                case 3 % z
                    % Draw z-grid
                    scolorname=searchcolor(id,get(ax,'ZColor'));
                    if strcmp(get(ax,'ZGrid'),'on')
                        if axlimz(1)~=axlimz(2)
                            xg_line_start = interp1([axlimz(1) axlimz(2)],[x(corners(k,2,p)) x(corners(k,3,p))],axztick);
                            yg_line_start = interp1([axlimz(1) axlimz(2)],[y(corners(k,2,p)) y(corners(k,3,p))],axztick);
                            xg_line_end = interp1([axlimz(1) axlimz(2)],[x(corners(k,4,p)) x(corners(k,5,p))],axztick);
                            yg_line_end = interp1([axlimz(1) axlimz(2)],[y(corners(k,4,p)) y(corners(k,5,p))],axztick);
                            for i = axzindex_inner
                                line2svg(fid,grouplabel,axpos,[xg_line_start(i) xg_line_end(i)],[yg_line_start(i) yg_line_end(i)],scolorname,gridlinestyle,linewidth);
                            end
                            if strcmp(get(ax,'ZTickMode'),'auto') && strcmp(get(ax,'ZMinorGrid'),'on') && ~isempty(minor_axztick)
                                xg_line_start = interp1([axlimz(1) axlimz(2)],[x(corners(k,2,p)) x(corners(k,3,p))],minor_axztick);
                                yg_line_start = interp1([axlimz(1) axlimz(2)],[y(corners(k,2,p)) y(corners(k,3,p))],minor_axztick);
                                xg_line_end = interp1([axlimz(1) axlimz(2)],[x(corners(k,4,p)) x(corners(k,5,p))],minor_axztick);
                                yg_line_end = interp1([axlimz(1) axlimz(2)],[y(corners(k,4,p)) y(corners(k,5,p))],minor_axztick);
                                for i = 1:length(xg_line_start)
                                    line2svg(fid,grouplabel,axpos,[xg_line_start(i) xg_line_end(i)],[yg_line_start(i) yg_line_end(i)],scolorname,minor_gridlinestyle,linewidth)
                                end  
                            end
                        end
                    end
                    if projection.xyplane == false
                        if strcmp(get(ax,'Box'),'on')
                            line2svg(fid,grouplabel,axpos,[x(corners(k,2,p)) x(corners(k,3,p))],[y(corners(k,2,p)) y(corners(k,3,p))],scolorname,'-',linewidth);
                            line2svg(fid,grouplabel,axpos,[x(corners(k,4,p)) x(corners(k,5,p))],[y(corners(k,4,p)) y(corners(k,5,p))],scolorname,'-',linewidth);
                        else
                            if strcmp(get(ax,'ZGrid'),'on')
                                line2svg(fid,grouplabel,axpos,[x(corners(k,2,p)) x(corners(k,3,p))],[y(corners(k,2,p)) y(corners(k,3,p))],scolorname,gridlinestyle,linewidth);
                                line2svg(fid,grouplabel,axpos,[x(corners(k,4,p)) x(corners(k,5,p))],[y(corners(k,4,p)) y(corners(k,5,p))],scolorname,gridlinestyle,linewidth);
                            end
                        end
                    end
            end
        end
    end
end
fprintf(fid,'    <g>\n');
axchild=get(ax,'Children');
group = axchild2svg(fid,id,axIdString,ax,group,paperpos,axchild,axpos,groupax,projection);
fprintf(fid,'    </g>\n');
if strcmp(get(ax,'Visible'),'on')
    fprintf(fid,'    <g>\n');
    % Search axis for labeling
    if projection.xyplane
        [x_axis_point_value, x_axis_point_index_top] = min(y);
        [x_axis_point_value, x_axis_point_index_bottom] = max(y);
        if strcmp(get(ax,'Box'),'on')
            if strcmp(get(ax,'XAxisLocation'),'top')
                x_axis_point_index = [x_axis_point_index_top x_axis_point_index_bottom];
            else
                x_axis_point_index = [x_axis_point_index_bottom x_axis_point_index_top];
            end
        else
            if strcmp(get(ax,'XAxisLocation'),'top')
                x_axis_point_index = x_axis_point_index_top;
            else
                x_axis_point_index = x_axis_point_index_bottom;
            end
        end
        [y_axis_point_value, y_axis_point_index_left] = min(x);
        [y_axis_point_value, y_axis_point_index_right] = max(x);
        if strcmp(get(ax,'Box'),'on')
            if strcmp(get(ax,'YAxisLocation'),'right')
                y_axis_point_index = [y_axis_point_index_right y_axis_point_index_left];
            else
                y_axis_point_index = [y_axis_point_index_left y_axis_point_index_right];
            end
        else
            if strcmp(get(ax,'YAxisLocation'),'right')
                y_axis_point_index = y_axis_point_index_right;
            else
                y_axis_point_index = y_axis_point_index_left;
            end
        end
        [z_axis_point_value, z_axis_point_index] = min(x);   
    else
        [x_axis_point_value, x_axis_point_index] = max(y);
        [y_axis_value, y_axis_point_index] = max(y);
        [z_axis_point_value, z_axis_point_index] = min(x);   
    end
    scolorname=searchcolor(id,get(ax,'XColor'));
    % Draw 'box' of x-axis
    if projection.xyplane == false
        if strcmp(get(ax,'Box'),'on')
            edge_line_index = [edge_opposite(most_back_edge_index) edge_neighbours(edge_opposite(most_back_edge_index),1)];
            line2svg(fid,grouplabel,axpos,x(edge_line_index),y(edge_line_index),scolorname,'-',linewidth)
        end
    end
    % Draw x-tick labels
    if (strcmp(get(ax,'XTickLabelMode'),'auto') && strcmp(get(ax,'XScale'),'log'))
        exponent = 1;
    else
        exponent = 0;
    end
    % Draw x-tick marks
    if (ticklength(1) ~= 0)
        if axlimx(1)~=axlimx(2)
            if (nomx(x_axis_point_index(1)))
                lim = [axlimx(2) axlimx(1)];    
            else
                lim = [axlimx(1) axlimx(2)];
            end
            x_label_end1 = interp1([0 1],[x(x_axis_point_index(1)) x(edge_neighbours(x_axis_point_index(1),2))],label_distance,'linear','extrap');
            y_label_end1 = interp1([0 1],[y(x_axis_point_index(1)) y(edge_neighbours(x_axis_point_index(1),2))],label_distance,'linear','extrap');
            x_label_end2 = interp1([0 1],[x(edge_neighbours(x_axis_point_index(1),1)) x(edge_neighbours(edge_neighbours(x_axis_point_index(1),1),2))],label_distance,'linear','extrap');
            y_label_end2 = interp1([0 1],[y(edge_neighbours(x_axis_point_index(1),1)) y(edge_neighbours(edge_neighbours(x_axis_point_index(1),1),2))],label_distance,'linear','extrap');
            xg_label_end = interp1(lim,[x_label_end1 x_label_end2],axxtick);
            yg_label_end = interp1(lim,[y_label_end1 y_label_end2],axxtick);            
            for k = 1:length(x_axis_point_index)
                x_tick_end1 = interp1([0 1],[x(x_axis_point_index(k)) x(edge_neighbours(x_axis_point_index(k),2))],ticklength*tick_ratio(1),'linear','extrap');
                y_tick_end1 = interp1([0 1],[y(x_axis_point_index(k)) y(edge_neighbours(x_axis_point_index(k),2))],ticklength*tick_ratio(1),'linear','extrap');
                x_tick_end2 = interp1([0 1],[x(edge_neighbours(x_axis_point_index(k),1)) x(edge_neighbours(edge_neighbours(x_axis_point_index(k),1),2))],ticklength*tick_ratio(1),'linear','extrap');
                y_tick_end2 = interp1([0 1],[y(edge_neighbours(x_axis_point_index(k),1)) y(edge_neighbours(edge_neighbours(x_axis_point_index(k),1),2))],ticklength*tick_ratio(1),'linear','extrap');
                xg_line_start = interp1(lim,[x(x_axis_point_index(k)) x(edge_neighbours(x_axis_point_index(k),1))],axxtick);
                yg_line_start = interp1(lim,[y(x_axis_point_index(k)) y(edge_neighbours(x_axis_point_index(k),1))],axxtick);
                xg_line_end = interp1(lim,[x_tick_end1 x_tick_end2],axxtick);
                yg_line_end = interp1(lim,[y_tick_end1 y_tick_end2],axxtick);
                for i = valid_xsticks
                    line2svg(fid,grouplabel,axpos,[xg_line_start(i) xg_line_end(i)],[yg_line_start(i) yg_line_end(i)],scolorname,'-',linewidth)
                end
                line2svg(fid,grouplabel,axpos,[x(x_axis_point_index(k)) x(edge_neighbours(x_axis_point_index(k),1))],[y(x_axis_point_index(k)) y(edge_neighbours(x_axis_point_index(k),1))],scolorname,'-',linewidth)
            end
            if ~isempty(axlabelx)
                % Note: 3D plot do not support the property XAxisLocation
                % setting 'top'.
                if strcmp(get(ax,'XAxisLocation'),'top') && (projection.xyplane == true)
                    for i = axxindex
                        label2svg(fid,grouplabel,axpos,ax,xg_label_end(i),yg_label_end(i),convertString(axlabelx(i,:)),'Center',0,'bottom',1,paperpos,scolorname,exponent);
                    end
                else
                    for i = axxindex
                        label2svg(fid,grouplabel,axpos,ax,xg_label_end(i),yg_label_end(i),convertString(axlabelx(i,:)),'Center',0,'top',1,paperpos,scolorname,exponent);
                    end
                end
            end
        end
    end
    scolorname=searchcolor(id,get(ax,'YColor'));
    % Draw 'box' of y-axis
    if projection.xyplane == false
        if strcmp(get(ax,'Box'),'on')
            edge_line_index = [edge_opposite(most_back_edge_index) edge_neighbours(edge_opposite(most_back_edge_index),2)];
            line2svg(fid,grouplabel,axpos,x(edge_line_index),y(edge_line_index),scolorname,'-',linewidth)
        end
    end
    % Draw y-tick labels
    if (strcmp(get(ax,'YTickLabelMode'),'auto') && strcmp(get(ax,'YScale'),'log'))
        exponent = 1;
    else
        exponent = 0;
    end
    % Draw y-tick marks
    if (ticklength(1) ~= 0)
        if axlimy(1)~=axlimy(2)
            if (nomx(y_axis_point_index(1)))
                lim = [axlimy(2) axlimy(1)];    
            else
                lim = [axlimy(1) axlimy(2)];
            end
            x_label_end1 = interp1([0 1],[x(y_axis_point_index(1)) x(edge_neighbours(y_axis_point_index(1),1))],label_distance,'linear','extrap');
            y_label_end1 = interp1([0 1],[y(y_axis_point_index(1)) y(edge_neighbours(y_axis_point_index(1),1))],label_distance,'linear','extrap');
            x_label_end2 = interp1([0 1],[x(edge_neighbours(y_axis_point_index(1),2)) x(edge_neighbours(edge_neighbours(y_axis_point_index(1),2),1))],label_distance,'linear','extrap');
            y_label_end2 = interp1([0 1],[y(edge_neighbours(y_axis_point_index(1),2)) y(edge_neighbours(edge_neighbours(y_axis_point_index(1),2),1))],label_distance,'linear','extrap');
            xg_label_end = interp1(lim,[x_label_end1 x_label_end2],axytick);
            yg_label_end = interp1(lim,[y_label_end1 y_label_end2],axytick);            
            for k = 1:length(y_axis_point_index)
                x_tick_end1 = interp1([0 1],[x(y_axis_point_index(k)) x(edge_neighbours(y_axis_point_index(k),1))],ticklength*tick_ratio(2),'linear','extrap');
                y_tick_end1 = interp1([0 1],[y(y_axis_point_index(k)) y(edge_neighbours(y_axis_point_index(k),1))],ticklength*tick_ratio(2),'linear','extrap');
                x_tick_end2 = interp1([0 1],[x(edge_neighbours(y_axis_point_index(k),2)) x(edge_neighbours(edge_neighbours(y_axis_point_index(k),2),1))],ticklength*tick_ratio(2),'linear','extrap');
                y_tick_end2 = interp1([0 1],[y(edge_neighbours(y_axis_point_index(k),2)) y(edge_neighbours(edge_neighbours(y_axis_point_index(k),2),1))],ticklength*tick_ratio(2),'linear','extrap');
                xg_line_start = interp1(lim,[x(y_axis_point_index(k)) x(edge_neighbours(y_axis_point_index(k),2))],axytick);
                yg_line_start = interp1(lim,[y(y_axis_point_index(k)) y(edge_neighbours(y_axis_point_index(k),2))],axytick);
                xg_line_end = interp1(lim,[x_tick_end1 x_tick_end2],axytick);
                yg_line_end = interp1(lim,[y_tick_end1 y_tick_end2],axytick);
                for i = valid_ysticks
                    line2svg(fid,grouplabel,axpos,[xg_line_start(i) xg_line_end(i)],[yg_line_start(i) yg_line_end(i)],scolorname,'-',linewidth)
                end
                line2svg(fid,grouplabel,axpos,[x(y_axis_point_index(k)) x(edge_neighbours(y_axis_point_index(k),2))],[y(y_axis_point_index(k)) y(edge_neighbours(y_axis_point_index(k),2))],scolorname,'-',linewidth)
            end
            if ~isempty(axlabely)
                % Note: 3D plot do not support the property YAxisLocation
                % setting 'right'.
                if (projection.xyplane == true)
                    if strcmp(get(ax,'YAxisLocation'),'right')
                        for i = axyindex
                            label2svg(fid,grouplabel,axpos,ax,xg_label_end(i),yg_label_end(i),convertString(axlabely(i,:)),'Left',0,'middle',1,paperpos,scolorname,exponent);
                        end
                   else
                        for i = axyindex
                            label2svg(fid,grouplabel,axpos,ax,xg_label_end(i),yg_label_end(i),convertString(axlabely(i,:)),'Right',0,'middle',1,paperpos,scolorname,exponent);
                        end
                    end
                else
                    for i = axyindex
                        label2svg(fid,grouplabel,axpos,ax,xg_label_end(i),yg_label_end(i),convertString(axlabely(i,:)),'Center',0,'top',1,paperpos,scolorname,exponent);
                    end
                end
            end
        end
    end
    scolorname=searchcolor(id,get(ax,'ZColor'));
    % Draw 'box' of z-axis
    if projection.xyplane == false
        if strcmp(get(ax,'Box'),'on')
            edge_line_index = [edge_opposite(most_back_edge_index) edge_neighbours(edge_opposite(most_back_edge_index),3)];
            line2svg(fid,grouplabel,axpos,x(edge_line_index),y(edge_line_index),scolorname,'-',linewidth)
        end
    end
    if (strcmp(get(ax,'ZTickLabelMode'),'auto') && strcmp(get(ax,'ZScale'),'log'))
        exponent = 1;
    else
        exponent = 0;
    end
    % Draw z-tick marks
    if (ticklength(1) ~= 0)
        if axlimz(1)~=axlimz(2)
            if (nomz(z_axis_point_index(1)))
                lim = [axlimz(2) axlimz(1)];    
            else
                lim = [axlimz(1) axlimz(2)];
            end
            x_tick_end1 = interp1([0 1],[x(z_axis_point_index) x(edge_neighbours(z_axis_point_index,2))],ticklength*tick_ratio(3),'linear','extrap');
            y_tick_end1 = interp1([0 1],[y(z_axis_point_index) y(edge_neighbours(z_axis_point_index,2))],ticklength*tick_ratio(3),'linear','extrap');
            x_tick_end2 = interp1([0 1],[x(edge_neighbours(z_axis_point_index,3)) x(edge_neighbours(edge_neighbours(z_axis_point_index,3),2))],ticklength*tick_ratio(3),'linear','extrap');
            y_tick_end2 = interp1([0 1],[y(edge_neighbours(z_axis_point_index,3)) y(edge_neighbours(edge_neighbours(z_axis_point_index,3),2))],ticklength*tick_ratio(3),'linear','extrap');
            x_label_end1 = interp1([0 1],[x(z_axis_point_index) x(edge_neighbours(z_axis_point_index,2))],label_distance,'linear','extrap');
            y_label_end1 = interp1([0 1],[y(z_axis_point_index) y(edge_neighbours(z_axis_point_index,2))],label_distance,'linear','extrap');
            x_label_end2 = interp1([0 1],[x(edge_neighbours(z_axis_point_index,3)) x(edge_neighbours(edge_neighbours(z_axis_point_index,3),2))],label_distance,'linear','extrap');
            y_label_end2 = interp1([0 1],[y(edge_neighbours(z_axis_point_index,3)) y(edge_neighbours(edge_neighbours(z_axis_point_index,3),2))],label_distance,'linear','extrap');
            xg_line_start = interp1(lim,[x(z_axis_point_index) x(edge_neighbours(z_axis_point_index,3))],axztick);
            yg_line_start = interp1(lim,[y(z_axis_point_index) y(edge_neighbours(z_axis_point_index,3))],axztick);
            xg_line_end = interp1(lim,[x_tick_end1 x_tick_end2],axztick);
            yg_line_end = interp1(lim,[y_tick_end1 y_tick_end2],axztick);
            xg_label_end = interp1(lim,[x_label_end1 x_label_end2],axztick);
            yg_label_end = interp1(lim,[y_label_end1 y_label_end2],axztick);            
            for i = valid_zsticks
                line2svg(fid,grouplabel,axpos,[xg_line_start(i) xg_line_end(i)],[yg_line_start(i) yg_line_end(i)],scolorname,'-',linewidth)
            end
            line2svg(fid,grouplabel,axpos,[x(z_axis_point_index) x(edge_neighbours(z_axis_point_index,3))],[y(z_axis_point_index) y(edge_neighbours(z_axis_point_index,3))],scolorname,'-',linewidth)
            if ~isempty(axlabelz)
                for i = axzindex
                    label2svg(fid,grouplabel,axpos,ax,xg_label_end(i),yg_label_end(i),convertString(axlabelz(i,:)),'Right',0,'middle',1,paperpos,scolorname,exponent);
                end
            end
        end
    end
    exponent2svg(fid,groupax,axpos,paperpos,ax)
    fprintf(fid,'    </g>\n');
end
fprintf(fid,'  </g>\n');
set(ax,'Units',originalAxesUnits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% take any axis children and create objects for them
function group=axchild2svg(fid,id,axIdString,ax,group,paperpos,axchild,axpos,groupax,projection)
global colorname
global PLOT2SVG_globals
for i=length(axchild):-1:1
    if strcmp(get(axchild(i), 'Visible'), 'off')
        % do nothing
    elseif strcmp(get(axchild(i),'Type'),'line')
        scolorname=searchcolor(id,get(axchild(i),'Color'));
        linestyle=get(axchild(i),'LineStyle');
        linewidth=get(axchild(i),'LineWidth');
        marker=get(axchild(i),'Marker');
        markeredgecolor=get(axchild(i),'MarkerEdgeColor');
        if ischar(markeredgecolor)
            switch markeredgecolor
                case 'none',markeredgecolorname='none';
                otherwise,markeredgecolorname=scolorname;  % if markeredgecolorname is 'auto' or something else set the markeredgecolorname to the line color
            end    
        else    
            markeredgecolorname=searchcolor(id,markeredgecolor);
        end
        markerfacecolor=get(axchild(i),'MarkerFaceColor');
        if ischar(markerfacecolor)
            switch markerfacecolor
                case 'none',markerfacecolorname='none';
                otherwise,markerfacecolorname=scolorname;  % if markerfacecolorname is 'auto' or something else set the markerfacecolorname to the line color
            end
        else
            markerfacecolorname=searchcolor(id,markerfacecolor);
        end
        markersize=get(axchild(i),'MarkerSize')/1.5;
        linex = get(axchild(i),'XData');
        linex = linex(:)'; % Octave stores the data in a column vector
        if strcmp(get(ax,'XScale'),'log')
            linex(find(linex<=0)) = NaN;
            linex=log10(linex);
        end
        liney=get(axchild(i),'YData');
        liney = liney(:)'; % Octave stores the data in a column vector        
        if strcmp(get(ax,'YScale'),'log')
            liney(find(liney<=0)) = NaN;
            liney=log10(liney);
        end
        linez=get(axchild(i),'ZData');
        linez = linez(:)'; % Octave stores the data in a column vector        
        if isempty(linez)
            linez = zeros(size(linex));    
        end
        if strcmp(get(ax,'ZScale'),'log')
            linez(find(linez<=0)) = NaN;
            linez=log10(linez);
        end
        % put a line into a group with its markers
        if strcmp(get(axchild(i),'Clipping'),'on')
            fprintf(fid,'<g clip-path="url(#%s)">\n',axIdString);
        else
            fprintf(fid,'<g>\n');
        end
        [x,y,z] = project(linex,liney,linez,projection);
        x = (x*axpos(3)+axpos(1))*paperpos(3);
        y = (1-(y*axpos(4)+axpos(2)))*paperpos(4);
        line2svg(fid,groupax,axpos,x,y,scolorname,linestyle,linewidth)
        % put the markers into a subgroup of the lines
        fprintf(fid,'<g>\n');
        switch marker
            case 'none';
            case '.',group=group+1;circle2svg(fid,group,axpos,x,y,markersize*0.25,'none',markeredgecolorname,linewidth);
            case 'o',group=group+1;circle2svg(fid,group,axpos,x,y,markersize*0.75,markeredgecolorname,markerfacecolorname,linewidth);
            case '+',group=group+1;patch2svg(fid,group,axpos,x'*ones(1,5)+ones(length(linex),1)*[-1 1 NaN 0 0]*markersize,y'*ones(1,5)+ones(length(liney),1)*[0 0 NaN -1 1]*markersize,markeredgecolorname,'-',linewidth,markeredgecolorname, 1, 1);   
            case '*',group=group+1;patch2svg(fid,group,axpos,x'*ones(1,11)+ones(length(linex),1)*[-1 1 NaN 0 0 NaN -0.7 0.7 NaN -0.7 0.7]*markersize,y'*ones(1,11)+ones(length(liney),1)*[0 0 NaN -1 1 NaN 0.7 -0.7 NaN -0.7 0.7]*markersize,markeredgecolorname,'-',linewidth,markeredgecolorname, 1, 1);
            case 'x',group=group+1;patch2svg(fid,group,axpos,x'*ones(1,5)+ones(length(linex),1)*[-0.7 0.7 NaN -0.7 0.7]*markersize,y'*ones(1,5)+ones(length(liney),1)*[0.7 -0.7 NaN -0.7 0.7]*markersize,markeredgecolorname,'-',linewidth,markeredgecolorname, 1, 1);
            %% Octave keeps s, d, p and h in the HandleGraphics object, for the square, diamond, pentagram, and hexagram markers, respectively -- Jakob Malm
            case {'square', 's'},group=group+1;patch2svg(fid,group,axpos,x'*ones(1,5)+ones(length(linex),1)*[-1 -1 1 1 -1]*markersize,y'*ones(1,5)+ones(length(liney),1)*[-1 1 1 -1 -1]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1);
            case {'diamond', 'd'},group=group+1;patch2svg(fid,group,axpos,x'*ones(1,5)+ones(length(linex),1)*[-0.7071 0 0.7071 0 -0.7071]*markersize,y'*ones(1,5)+ones(length(liney),1)*[0 1 0 -1 0]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1);
            case {'pentagram', 'p'},group=group+1;patch2svg(fid,group,axpos,...
                    x'*ones(1,11)+ones(length(linex),1)*[0 0.1180 0.5 0.1910 0.3090 0 -0.3090 -0.1910 -0.5 -0.1180 0]*1.3*markersize,...
                    y'*ones(1,11)+ones(length(liney),1)*[-0.5257 -0.1625 -0.1625 0.0621 0.4253 0.2008 0.4253 0.0621 -0.1625 -0.1625 -0.5257]*1.3*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1);
            case {'hexagram', 'h'},group=group+1;patch2svg(fid,group,axpos,...
                    x'*ones(1,13)+ones(length(linex),1)*[0 0.2309 0.6928 0.4619 0.6928 0.2309 0 -0.2309 -0.6928 -0.4619 -0.6928 -0.2309 0]*1*markersize,...
                    y'*ones(1,13)+ones(length(liney),1)*[0.8 0.4 0.4 0 -0.4 -0.4 -0.8 -0.4 -0.4 0 0.4 0.4 0.8]*1*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1);    
            case '^',group=group+1;patch2svg(fid,group,axpos,x'*ones(1,4)+ones(length(linex),1)*[-1 1 0 -1]*markersize,y'*ones(1,4)+ones(length(liney),1)*[0.577 0.577 -0.837 0.577]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1);
            case 'v',group=group+1;patch2svg(fid,group,axpos,x'*ones(1,4)+ones(length(linex),1)*[-1 1 0 -1]*markersize,y'*ones(1,4)+ones(length(liney),1)*[-0.577 -0.577 0.837 -0.577]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1);
            case '<',group=group+1;patch2svg(fid,group,axpos,x'*ones(1,4)+ones(length(linex),1)*[0.577 0.577 -0.837 0.577]*markersize,y'*ones(1,4)+ones(length(liney),1)*[-1 1 0 -1]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1);
            case '>',group=group+1;patch2svg(fid,group,axpos,x'*ones(1,4)+ones(length(linex),1)*[-0.577 -0.577 0.837 -0.577]*markersize,y'*ones(1,4)+ones(length(liney),1)*[-1 1 0 -1]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1);
        end
        % close the marker group
        fprintf(fid,'</g>\n');
        % close the line group
        fprintf(fid,'</g>\n');
    elseif strcmp(get(axchild(i),'Type'),'patch')
        flat_shading = 1;
        cmap=get(id,'Colormap');
        pointc=get(axchild(i),'FaceVertexCData');
        %pointc=get(axchild(i),'CData');
        % Scale color if scaled color mapping is turned on
        if strcmp(get(axchild(i),'CDataMapping'),'scaled')
            clim=get(ax,'CLim');
            pointc=(pointc-clim(1))/(clim(2)-clim(1))*(size(cmap,1)-1)+1;
        end
        % Limit index to smallest or biggest color index
        pointc=max(pointc,1);
        pointc=min(pointc,size(cmap,1));
        if ~ischar(get(axchild(i),'FaceAlpha'))
            face_opacity = get(axchild(i),'FaceAlpha');
        else
            face_opacity = 1.0;
        end
        if ~ischar(get(axchild(i),'EdgeAlpha'))
            edge_opacity = get(axchild(i),'EdgeAlpha');
        else
            edge_opacity = 1.0;
        end
        linestyle=get(axchild(i),'LineStyle');
        linewidth=get(axchild(i),'LineWidth');
        marker=get(axchild(i),'Marker');
        markeredgecolor=get(axchild(i),'MarkerEdgeColor');
        markersize=get(axchild(i),'MarkerSize')/1.5;
        points=get(axchild(i),'Vertices')';
        if strcmp(get(ax,'XScale'),'log')
            points(1,:)=log10(points(1,:));
        end
        if strcmp(get(ax,'YScale'),'log')
            points(2,:)=log10(points(2,:));
        end
        % TODO LogZ
        if size(points,1)==2
            [x,y,z] = project(points(1,:),points(2,:),zeros(size(points(1,:))),projection);    
        else
            [x,y,z] = project(points(1,:),points(2,:),points(3,:),projection);
        end
        x = (x*axpos(3)+axpos(1))*paperpos(3);
        y = (1-(y*axpos(4)+axpos(2)))*paperpos(4);
        faces = get(axchild(i),'Faces');
        face_index = 1:size(faces,1);
        if size(points,1)==3;
            [z,face_index]=sort(sum(z(faces(:,:)),2));
            faces=faces(face_index,:);
        end
        if strcmp(get(axchild(i),'Clipping'),'on')
            fprintf(fid,'<g clip-path="url(#%s)">\n',axIdString);
        else
            fprintf(fid,'<g>\n');
        end
        for p=1:size(faces,1)
            if ischar(get(axchild(i),'FaceColor'))
                if strcmp(get(axchild(i),'FaceColor'),'texturemap')
                    facecolorname='none';   % TO DO: texture map
                elseif strcmp(get(axchild(i),'FaceColor'),'none')
                    facecolorname='none';
                else
                    if size(pointc,1)==1
                        facecolor = pointc;    
                    elseif size(pointc,1)==size(faces,1)
                        if strcmp(get(axchild(i),'FaceColor'),'flat')
                            facecolor = pointc(face_index(p),:);
                        else
                            facecolor = pointc(face_index(p),:);
                            cdata = pointc(face_index(p),:);     % TO DO: color interpolation
                            flat_shading = 0;
                        end
                    elseif size(pointc,1)==size(points,2)
                        if strcmp(get(axchild(i),'FaceColor'),'flat')
                            facecolor = pointc(faces(p,1),:);
                        else
                            facecolor = pointc(faces(p,1));
                            cdata = pointc(faces(p,:),:);
                            flat_shading = 0;
                        end
                    else
                        error('Unsupported color handling for patches.');    
                    end
                    if ~isnan(facecolor)
                        if size(facecolor,2)==1
                            facecolorname = ['#' colorname(ceil(facecolor),:)];
                        else
                            if strcmp(get(axchild(i),'FaceColor'),'flat')  % Bugfix 27.01.2008
                                facecolorname = searchcolor(id,facecolor/64);
                            else
                                facecolorname = searchcolor(id,facecolor);    
                            end
                        end
                    else
                        facecolorname='none';
                    end
                end
            else
                facecolorname = searchcolor(id,get(axchild(i),'FaceColor'));       
            end
            if ischar(get(axchild(i),'EdgeColor'))
                if strcmp(get(axchild(i),'EdgeColor'),'none')
                    edgecolorname = 'none';
                else
                    if size(pointc,1)==1
                        edgecolor = pointc;    
                    elseif size(pointc,1)==size(faces,1)
                        edgecolor = pointc(p,:);
                    elseif size(pointc,1)==size(points,2)
                        if strcmp(get(axchild(i),'EdgeColor'),'flat')
                            edgecolor = pointc(faces(p,1));
                        else
                            edgecolor = pointc(faces(p,1));     % TO DO: color interpolation
                        end
                    else
                        error('Unsupported color handling for patches.');    
                    end
                    if ~isnan(edgecolor)
                        if size(edgecolor,2)==1
                            edgecolorname = ['#' colorname(ceil(edgecolor),:)];
                        else
                            if strcmp(get(axchild(i),'EdgeColor'),'flat')   % Bugfix 27.01.2008
                                edgecolorname = searchcolor(id,edgecolor/64);
                            else
                                edgecolorname = searchcolor(id,edgecolor);    
                            end
                        end
                    else
                        edgecolorname = 'none';
                    end
                end
            else
                edgecolorname = searchcolor(id,get(axchild(i),'EdgeColor'));       
            end
            if flat_shading
                patch2svg(fid, group, axpos, x(faces(p,:)), y(faces(p,:)), facecolorname, linestyle, linewidth, edgecolorname, face_opacity, edge_opacity)
            else
                %patch2svg(fid, group, axpos, x(faces(p,:)), y(faces(p,:)), facecolorname, linestyle, linewidth, edgecolorname, face_opacity, edge_opacity)
                gouraud_patch2svg(fid, group, axpos, x(faces(p,:)), y(faces(p,:)), cdata, linestyle, linewidth, edgecolorname, face_opacity, edge_opacity,id)
            end
            if ~strcmp(marker, 'none')
                xmarker = x(faces(p,:));
                ymarker = y(faces(p,:));
                % put the markers into a subgroup of the lines
                fprintf(fid,'<g>\n');
                if ischar(markeredgecolor)
                    switch markeredgecolor
                        case 'none',markeredgecolorname='none';
                        otherwise,markeredgecolorname=edgecolorname;  % if markeredgecolorname is 'auto' or something else set the markeredgecolorname to the line color
                    end    
                else    
                    markeredgecolorname=searchcolor(id,markeredgecolor);
                end
                markerfacecolor=get(axchild(i),'MarkerFaceColor');
                if ischar(markerfacecolor)
                    switch markerfacecolor
                        case 'none',markerfacecolorname='none';
                        otherwise,markerfacecolorname=edgecolorname;  % if markerfacecolorname is 'auto' or something else set the markerfacecolorname to the line color
                    end
                else
                    markerfacecolorname=searchcolor(id,markerfacecolor);
                end
                
                switch marker
                    case 'none';
                    case '.',group=group+1;circle2svg(fid,group,axpos,xmarker,ymarker,markersize*0.25,'none',markeredgecolorname,linewidth);
                    case 'o',group=group+1;circle2svg(fid,group,axpos,xmarker,ymarker,markersize*0.75,markeredgecolorname,markerfacecolorname,linewidth);
                    case '+',group=group+1;patch2svg(fid,group,axpos,xmarker'*ones(1,5)+ones(length(linex),1)*[-1 1 NaN 0 0]*markersize,ymarker'*ones(1,5)+ones(length(liney),1)*[0 0 NaN -1 1]*markersize,markeredgecolorname,'-',linewidth,markeredgecolorname, 1, 1);   
                    case '*',group=group+1;patch2svg(fid,group,axpos,xmarker'*ones(1,11)+ones(length(linex),1)*[-1 1 NaN 0 0 NaN -0.7 0.7 NaN -0.7 0.7]*markersize,ymarker'*ones(1,11)+ones(length(liney),1)*[0 0 NaN -1 1 NaN 0.7 -0.7 NaN -0.7 0.7]*markersize,markeredgecolorname,'-',linewidth,markeredgecolorname, 1, 1);
                    case 'x',group=group+1;patch2svg(fid,group,axpos,xmarker'*ones(1,5)+ones(length(linex),1)*[-0.7 0.7 NaN -0.7 0.7]*markersize,ymarker'*ones(1,5)+ones(length(liney),1)*[0.7 -0.7 NaN -0.7 0.7]*markersize,markeredgecolorname,'-',linewidth,markeredgecolorname, 1, 1);
                    %% Octave keeps s, d, p and h in the HandleGraphics object, for the square, diamond, pentagram, and hexagram markers, respectively -- Jakob Malm    
                    case {'square', 's'},group=group+1;patch2svg(fid,group,axpos,xmarker'*ones(1,5)+ones(length(linex),1)*[-1 -1 1 1 -1]*markersize,ymarker'*ones(1,5)+ones(length(liney),1)*[-1 1 1 -1 -1]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1);
                    case {'diamond', 'd'},group=group+1;patch2svg(fid,group,axpos,xmarker'*ones(1,5)+ones(length(linex),1)*[-0.7071 0 0.7071 0 -0.7071]*markersize,ymarker'*ones(1,5)+ones(length(liney),1)*[0 1 0 -1 0]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1);
                    case {'pentagram', 'p'},group=group+1;patch2svg(fid,group,axpos,...
                            xmarker'*ones(1,11)+ones(length(linex),1)*[0 0.1180 0.5 0.1910 0.3090 0 -0.3090 -0.1910 -0.5 -0.1180 0]*1.3*markersize,...
                            ymarker'*ones(1,11)+ones(length(liney),1)*[-0.5257 -0.1625 -0.1625 0.0621 0.4253 0.2008 0.4253 0.0621 -0.1625 -0.1625 -0.5257]*1.3*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1);
                    case {'hexagram', 'h'},group=group+1;patch2svg(fid,group,axpos,...
                            xmarker'*ones(1,13)+ones(length(linex),1)*[0 0.2309 0.6928 0.4619 0.6928 0.2309 0 -0.2309 -0.6928 -0.4619 -0.6928 -0.2309 0]*1*markersize,...
                            ymarker'*ones(1,13)+ones(length(liney),1)*[0.8 0.4 0.4 0 -0.4 -0.4 -0.8 -0.4 -0.4 0 0.4 0.4 0.8]*1*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1);    
                    case '^',group=group+1;patch2svg(fid,group,axpos,xmarker'*ones(1,4)+ones(length(linex),1)*[-1 1 0 -1]*markersize,ymarker'*ones(1,4)+ones(length(liney),1)*[0.577 0.577 -0.837 0.577]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1);
                    case 'v',group=group+1;patch2svg(fid,group,axpos,xmarker'*ones(1,4)+ones(length(linex),1)*[-1 1 0 -1]*markersize,ymarker'*ones(1,4)+ones(length(liney),1)*[-0.577 -0.577 0.837 -0.577]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1);
                    case '<',group=group+1;patch2svg(fid,group,axpos,xmarker'*ones(1,4)+ones(length(linex),1)*[0.577 0.577 -0.837 0.577]*markersize,ymarker'*ones(1,4)+ones(length(liney),1)*[-1 1 0 -1]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1);
                    case '>',group=group+1;patch2svg(fid,group,axpos,xmarker'*ones(1,4)+ones(length(linex),1)*[-0.577 -0.577 0.837 -0.577]*markersize,ymarker'*ones(1,4)+ones(length(liney),1)*[-1 1 0 -1]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1);
                end
                % close the marker group
                fprintf(fid,'</g>\n');
            end
        end
        fprintf(fid,'</g>\n');
    elseif strcmp(get(axchild(i),'Type'),'surface')
        flat_shading = 1;
        cmap=get(id,'Colormap');
        [faces,points,pointc,alpha]=surface2patch(axchild(i));
        points=points';
        % Scale color if scaled color mapping is turned on
        if strcmp(get(axchild(i),'CDataMapping'),'scaled')
            clim=get(ax,'CLim');
            pointc=(pointc-clim(1))/(clim(2)-clim(1))*(size(cmap,1)-1)+1;
        end
        % Limit index to smallest or biggest color index
        pointc=max(pointc,1);
        if ~ischar(get(axchild(i),'FaceAlpha'))
            face_opacity = get(axchild(i),'FaceAlpha');
        elseif strcmp(get(axchild(i),'FaceAlpha'),'flat')
            face_opacity = alpha;
            switch get(axchild(i),'AlphaDataMapping')
                case {'direct'}
                    face_opacity = 1.0; % TODO
                case {'scaled'}
                    alim=get(ax,'ALim');
                    face_opacity=(face_opacity-alim(1))/(alim(2)-alim(1));
                case {'none'}
                    % Clip alpha data
                    face_opacity = min(1, face_opacity);
                    face_opacity = max(0, face_opacity);
                otherwise
                    error(['Unsupported AlphaDataMapping identifier ''' get(axchild(i),'AlphaDataMapping') '''.']);
            end
        else
            face_opacity = 1.0;
        end
        if ~ischar(get(axchild(i),'EdgeAlpha'))
            edge_opacity = get(axchild(i),'EdgeAlpha');
        else
            edge_opacity = 1.0;
        end
        pointc=min(pointc,size(cmap,1));
        linestyle=get(axchild(i),'LineStyle');
        linewidth=get(axchild(i),'LineWidth');
        if strcmp(get(ax,'XScale'),'log')
            points(1,:)=log10(points(1,:));
        end
        if strcmp(get(ax,'YScale'),'log')
            points(2,:)=log10(points(2,:));
        end
        if size(points,1)==3
            if strcmp(get(ax,'ZScale'),'log')
                points(3,:)=log10(points(3,:));
            end   
        end
        if size(points,1)==3
            [x,y,z] = project(points(1,:),points(2,:),points(3,:),projection); 
        else
            [x,y,z] = project(points(1,:),points(2,:),zeros(size(points(1,:))),projection); 
        end
        x = (x*axpos(3)+axpos(1))*paperpos(3);
        y = (1-(y*axpos(4)+axpos(2)))*paperpos(4);
        face_index = 1:size(faces,1);
        if size(points,1)==3;
            [z,face_index]=sort(sum(z(faces(:,:)),2));
            faces=faces(face_index,:);
        end
        if strcmp(get(axchild(i),'Clipping'),'on')
            fprintf(fid,'<g clip-path="url(#%s)">\n',axIdString);
        else
            fprintf(fid,'<g>\n');
        end
        for p=1:size(faces,1)
            if ischar(get(axchild(i),'FaceColor'))
                if strcmp(get(axchild(i),'FaceColor'),'texturemap')
                    facecolorname='none';   % TO DO: texture map
                elseif strcmp(get(axchild(i),'FaceColor'),'none')
                    facecolorname='none';
                else
                    if size(pointc,1)==1
                        facecolor = pointc;    
                    elseif size(pointc,1)==size(faces,1)
                        facecolor = pointc(face_index(p),:);
                    elseif size(pointc,1)==size(points,2)
                        if strcmp(get(axchild(i),'FaceColor'),'flat')
                            facecolor = pointc(faces(p,1));
                        else
                            facecolor = pointc(faces(p,1));
                            cdata = pointc(faces(p,:));
                            flat_shading = 0;
                        end
                    else
                        error('Unsupported color handling for patches.');    
                    end
                    if ~isnan(facecolor)
                        if size(facecolor,2)==1
                            facecolorname = ['#' colorname(ceil(facecolor),:)];
                        else
                            facecolorname = searchcolor(id,facecolor);    
                        end
                    else
                        facecolorname='none';
                    end
                end
            else
                facecolorname = searchcolor(id,get(axchild(i),'FaceColor'));       
            end
            if size(face_opacity,1)==1
                face_opacity_value = face_opacity;
            elseif size(face_opacity,1)==size(faces,1)
                face_opacity_value = face_opacity(p,:);
            elseif size(face_opacity,1)==size(points,2)
                face_opacity_value = face_opacity(faces(p,1));
            else
                error('Unsupported face alpha value handling for patches.');
            end
            if ischar(get(axchild(i),'EdgeColor'))
                if strcmp(get(axchild(i),'EdgeColor'),'none')
                    edgecolorname = 'none';
                else
                    if size(pointc,1)==1
                        edgecolor = pointc;    
                    elseif size(pointc,1)==size(faces,1)
                        edgecolor = pointc(p,:);
                    elseif size(pointc,1)==size(points,2)
                        if strcmp(get(axchild(i),'EdgeColor'),'flat')
                            edgecolor = pointc(faces(p,1));
                        else
                            edgecolor = pointc(faces(p,1));     % TO DO: color interpolation
                        end
                    else
                        error('Unsupported color handling for patches.');    
                    end
                    if ~isnan(edgecolor)
                        if size(edgecolor,2)==1
                            edgecolorname = ['#' colorname(ceil(edgecolor),:)];
                        else
                            edgecolorname = searchcolor(id,edgecolor);    
                        end
                    else
                        edgecolorname = 'none';
                    end
                end
            else
                edgecolorname = searchcolor(id,get(axchild(i),'EdgeColor'));       
            end
            if flat_shading
                patch2svg(fid, group, axpos, x(faces(p,:)), y(faces(p,:)), facecolorname, linestyle, linewidth, edgecolorname,  face_opacity_value, edge_opacity)
            else
                gouraud_patch2svg(fid, group, axpos, x(faces(p,:)), y(faces(p,:)), cdata, linestyle, linewidth, edgecolorname, face_opacity_value, edge_opacity,id)
            end
        end
        fprintf(fid,'</g>\n');
    elseif strcmp(get(axchild(i),'Type'),'text')
        if strcmp(get(axchild(i),'Clipping'),'on')
            fprintf(fid,'<g clip-path="url(#%s)">\n',axIdString);
            text2svg(fid,1,axpos,paperpos,axchild(i),ax,projection)
            fprintf(fid,'</g>\n');
        else
            text2svg(fid,1,axpos,paperpos,axchild(i),ax,projection)
        end
    elseif strcmp(get(axchild(i),'Type'),'image')
        cmap=get(id,'Colormap');
        pointx=get(axchild(i),'XData');
        pointy=get(axchild(i),'YData');
        % If the XData is a vector we only use start and stop for the image
        if (size(pointx, 1) > 2) || (size(pointx, 2) > 1)
            pointx = [pointx(1) pointx(end)];   
        end
        if (size(pointy, 1) > 2) || (size(pointy, 2) > 1) 
            pointy = [pointy(1) pointy(end)];   
        end
        if (size(pointx, 1) > 1) && (size(pointy, 1) > 1)
            [x,y,z] = project(pointx,zeros(size(pointx)),zeros(size(pointx)),projection);
        else
            [x,y_dummy,z] = project(pointx,zeros(size(pointx)),zeros(size(pointx)),projection);
            [x_dummy,y,z] = project(zeros(size(pointy)),pointy,zeros(size(pointy)),projection);
        end
        pointc=get(axchild(i),'CData');
        %pointcclass = class(pointc);  % Bugfix proposed by Tom
        if strcmp(get(axchild(i),'CDataMapping'),'scaled')
            clim=get(ax,'CLim');
            pointc=(pointc-clim(1))/(clim(2)-clim(1))*(size(cmap,1) - 1) + 1; % Bugfix proposed by Tom
            %pointcclass = 'double'; % since range is now [0->size(cmap,1)-1]  % Bugfix proposed by Tom
        end
        data_aspect_ratio = get(ax,'DataAspectRatio');
        if length(x) == 2
            if size(pointc, 2) == 1
                halfwidthx = abs(x(2) - x(1)) * data_aspect_ratio(1);
            else
                halfwidthx = abs(x(2) - x(1))/(size(pointc,2) - 1);   
            end
        else
            halfwidthx = data_aspect_ratio(1);
        end
        if length(y) == 2
            if size(pointc, 1) == 1
                halfwidthy = abs(y(2) - y(1)) * data_aspect_ratio(2);
            else
                halfwidthy = abs(y(2) - y(1))/(size(pointc,1) - 1);   
            end
        else
            halfwidthy = data_aspect_ratio(2);
        end
        if length(pointx) > 1
            if xor(strcmp(get(ax,'XDir'),'reverse'), pointx(1) > pointx(2))
                if ndims(pointc) < 3
                    pointc=fliplr(pointc);
                elseif ndims(pointc) == 3
                    for j = [1:size(pointc,3)]
                        pointc(:,:,j)=fliplr(pointc(:,:,j));
                    end
                else
                    error('Invalid number of dimensions of data.');
                end
            end
        end
        if length(pointy) > 1
            if xor(strcmp(get(ax,'YDir'),'reverse'), pointy(1) > pointy(2))
                if ndims(pointc) < 3
                    pointc=flipud(pointc);
                elseif ndims(pointc) == 3
                    for j = [1:size(pointc,3)]
                        pointc(:,:,j)=flipud(pointc(:,:,j));
                    end
                else
                    error('Invalid number of dimensions of data.');
                end
            end
        end
        % pointc = cast(pointc,pointcclass);  % Bugfix proposed by Tom
        % Function 'cast' is not supported by old Matlab versions
        if (~isa(pointc, 'double') && ~isa(pointc, 'single'))
            if strcmp(get(axchild(i),'CDataMapping'),'scaled')
                pointc = double(pointc);
            else
                pointc = double(pointc) + 1;
            end
        end
        if ndims(pointc) ~= 3
            pointc = max(min(round(double(pointc)),size(cmap,1)),1);
        end
        CameraUpVector=get(ax,'CameraUpVector');
        filename = [PLOT2SVG_globals.basefilename sprintf('%03d',PLOT2SVG_globals.figurenumber) '.' PLOT2SVG_globals.pixelfiletype];
        PLOT2SVG_globals.figurenumber = PLOT2SVG_globals.figurenumber + 1;
        if isempty(PLOT2SVG_globals.basefilepath)
            current_path = pwd;
        else
            current_path = PLOT2SVG_globals.basefilepath;
        end
        if exist(fullfile(current_path,filename),'file')
            lastwarn('');
            delete(filename);
            if strcmp(lastwarn,'File not found or permission denied.')
                error('Cannot write image file. Make sure that no image is opened in an other program.')    
            end
        end
        if ndims(pointc) < 3
            pointc = flipud(pointc);
        elseif ndims(pointc) == 3
            for j = size(pointc,3)
                pointc(:,:,j)=flipud(pointc(:,:,j));
            end
        else
            error('Invalid number of dimensions of data.');
        end
        if ndims(pointc) == 3
            %pointc is not indexed
            if PLOT2SVG_globals.octave
                imwrite(fullfile(PLOT2SVG_globals.basefilepath,filename),pointc);
            else
                imwrite(pointc,fullfile(PLOT2SVG_globals.basefilepath,filename),PLOT2SVG_globals.pixelfiletype);
            end
        else
            %pointc is probably indexed
            if PLOT2SVG_globals.octave
                imwrite(fullfile(PLOT2SVG_globals.basefilepath,filename),pointc,cmap);
            else
                imwrite(pointc,cmap,fullfile(PLOT2SVG_globals.basefilepath,filename),PLOT2SVG_globals.pixelfiletype);
            end
        end
            lx=(size(pointc,2)*halfwidthx)*axpos(3)*paperpos(3);
        	ly=(size(pointc,1)*halfwidthy)*axpos(4)*paperpos(4);
        if strcmp(get(ax,'DataAspectRatioMode'),'manual')
            pointsx=((min(x) - halfwidthx/2)*axpos(3)+axpos(1))*paperpos(3);
            pointsy=(1-((max(y) + halfwidthy/2)*axpos(4)+axpos(2)))*paperpos(4);
        else
            pointsx=axpos(1)*paperpos(3);
            pointsy=(1-(axpos(4)+axpos(2)))*paperpos(4);
        end
        if strcmp(get(axchild(i),'Clipping'),'on')
            fprintf(fid,'<g clip-path="url(#%s)">\n',axIdString);
            fprintf(fid,'<image x="%0.3f" y="%0.3f" width="%0.3f" height="%0.3f" image-rendering="optimizeSpeed" preserveAspectRatio="none" xlink:href="%s" />\n', pointsx, pointsy, lx, ly, filename);
            fprintf(fid,'</g>\n');
        else
            fprintf(fid,'<image x="%0.3f" y="%0.3f" width="%0.3f" height="%0.3f" image-rendering="optimizeSpeed" preserveAspectRatio="none" xlink:href="%s" />\n', pointsx, pointsy, lx, ly, filename);
        end
    elseif strcmp(get(axchild(i),'Type'), 'hggroup')
        % handle group types (like error bars)
        % FIXME: they are not yet perfectly handled, there are more options
        % that are not used
        if strcmp(get(axchild(i),'Clipping'),'on')
            fprintf(fid,'<g clip-path="url(#%s)">\n',axIdString);
        else
            fprintf(fid, '<g>');
        end
        group=axchild2svg(fid,id,axIdString,ax,group,paperpos,get(axchild(i), 'Children'),axpos,groupax,projection);
        fprintf(fid, '</g>');
    elseif strcmp(get(axchild(i),'Type'), 'hgtransform')
        if strcmpi(get(axchild(i), 'Visible'), 'on')
            if strcmp(get(axchild(i),'Clipping'),'on')
                fprintf(fid,'<g clip-path="url(#%s)">\n',axIdString);
            else
                fprintf(fid, '<g>');
            end
            group=axchild2svg(fid,id,axIdString,ax,group,paperpos,get(axchild(i), 'Children'),axpos,groupax,projection);
            fprintf(fid, '</g>');
        end
    else
        disp(['   Warning: Unhandled child type: ' get(axchild(i),'Type')]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a patch (filled area)
function patch2svg(fid,group,axpos,xtot,ytot,scolorname,style,width, edgecolorname, face_opacity, edge_opacity)
for i=1:size(xtot,1)
    x=xtot(i,:);
    y=ytot(i,:);
    switch style
        case '--',pattern = 'stroke-dasharray="100pt,25pt"';
        case ':',pattern = 'stroke-dasharray="25pt,25pt"';
        case '-.',pattern = 'stroke-dasharray="100pt,25pt,25pt,25pt"';
        case 'none',pattern = 'stroke-dasharray="none"'; edge_opacity = 0.0;
        otherwise,pattern='stroke-dasharray="none"';   
    end 
    if (isnan(x)==zeros(size(x))&isnan(y)==zeros(size(y)))
        for j=1:20000:length(x)
            xx=x(j:min(length(x),j+19999));
            yy=y(j:min(length(y),j+19999));
            if ~strcmp(edgecolorname,'none') || ~strcmp(scolorname,'none')
                fprintf(fid,'      <polyline fill="%s" fill-opacity="%0.2f" stroke="%s" stroke-width="%0.1fpt" stroke-opacity="%0.2f" %s points="',...
                    scolorname, face_opacity, edgecolorname, width, edge_opacity, pattern);
                fprintf(fid,'%0.3f,%0.3f ',[xx;yy]);
                fprintf(fid,'"/>\n');        
            end
        end
    else
        parts=find(isnan(x)+isnan(y));
        if parts(1)~=1
            parts=[0 parts];
        end
        if parts(length(parts))~=length(x)
            parts=[parts length(x)+1];
        end
        for j=1:(length(parts)-1)
            xx=x((parts(j)+1):(parts(j+1)-1));
            yy=y((parts(j)+1):(parts(j+1)-1));
            if ~strcmp(edgecolorname,'none') || ~strcmp(scolorname,'none')
                if length(xx)~=0
                    fprintf(fid,'      <polyline fill="%s" fill-opacity="%0.2f" stroke="%s" stroke-width="%0.1fpt" stroke-opacity="%0.2f" %s points="',...
                        scolorname, face_opacity, edgecolorname, width, edge_opacity, pattern);
                    fprintf(fid,'%0.3f,%0.3f ',[xx;yy]);
                    fprintf(fid,'"/>\n');        
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a patch (filled area)
function gouraud_patch2svg(fid,group,axpos,xtot,ytot,cdata,style,width, edgecolorname, face_opacity, edge_opacity,id)
global colorname
for i=1:size(xtot,1)
    x=xtot(i,:);
    y=ytot(i,:);
    switch style
        case '--',pattern = 'stroke-dasharray="100pt,25pt"';
        case ':',pattern = 'stroke-dasharray="25pt,25pt"';
        case '-.',pattern = 'stroke-dasharray="100pt,25pt,25pt,25pt"';
        case 'none',pattern = 'stroke-dasharray="none"'; edge_opacity = 0.0;
        otherwise,pattern='stroke-dasharray="none"';   
    end 
    if (any(isnan(x)) || any(isnan(y)))
        fprintf('Warning: Found NaN in Gouraud patch.\n')
    else
        % If there are more than 2 edges always 3 edges are taken togehter
        % to form a triangle
        if length(x) > 2
            for j = 3:length(x)
                coord = [x([1 j-1 j]);y([1 j-1 j])];
                face_color = cdata(1,:);
                face_color2 = cdata(j-1,:);
                face_color3 = cdata(j,:);
                delta = coord(:,3)-coord(:,2);
                if det([delta (coord(:,1)-coord(:,2))]) ~= 0
                    if ~isnan(face_color)
                        IDstring1 = createId;
                        IDstring2 = createId;
                        if size(face_color2,2)==1
                            face_color_name2 = ['#' colorname(ceil(face_color2),:)];
                        else
                            face_color_name2 = searchcolor(id,face_color2/64);    
                        end
                        if size(face_color3,2)==1
                            face_color_name3 = ['#' colorname(ceil(face_color3),:)];
                        else
                            face_color_name3 = searchcolor(id,face_color3/64);    
                        end
                        grad_end=(delta)*(delta'*(coord(:,1)-coord(:,2)))/(delta'*delta) + coord(:,2);
                        if size(face_color,2)==1
                            face_color_name = ['#' colorname(ceil(face_color),:)];
                        else
                            face_color_name = searchcolor(id,face_color/64);    
                        end
                        fprintf(fid,'<defs>\n');
                        fprintf(fid,'<linearGradient id="%s" gradientUnits="userSpaceOnUse" x1="%0.3f" y1="%0.3f" x2="%0.3f" y2="%0.3f">\n',...
                            IDstring1, coord(1,2), coord(2,2), coord(1,3), coord(2,3));
                        fprintf(fid,'<stop offset="0" stop-color="%s" stop-opacity="1"/>\n',face_color_name2);
                        fprintf(fid,'<stop offset="1" stop-color="%s" stop-opacity="1"/>\n',face_color_name3);
                        fprintf(fid,'</linearGradient>\n');
                        fprintf(fid,'<linearGradient id="%s" gradientUnits="userSpaceOnUse" x1="%0.3f" y1="%0.3f" x2="%0.3f" y2="%0.3f">\n',...
                            IDstring2, coord(1,1), coord(2,1), grad_end(1), grad_end(2));
                        fprintf(fid,'<stop offset="0" stop-color="%s" stop-opacity="1"/>\n',face_color_name);
                        fprintf(fid,'<stop offset="1" stop-color="%s" stop-opacity="0"/>\n',face_color_name);
                        fprintf(fid,'</linearGradient>\n');
                        fprintf(fid,'</defs>\n');    
                        % Open group
                        temp_string = sprintf('%0.3f,%0.3f ',coord);
                        fprintf(fid,'<g opacity="%0.2f">\n',face_opacity);
                        fprintf(fid,'<polyline fill="url(#%s)" stroke="none" points="%s"/>\n',IDstring1,temp_string);
                        fprintf(fid,'<polyline fill="url(#%s)" stroke="none" points="%s"/>\n',IDstring2,temp_string);
                        % Close group
                        fprintf(fid,'</g>\n');
                    end
                end
            end
        end
        % Last we draw the line around the patch
        if ~strcmp(edgecolorname,'none')
            fprintf(fid,'<polygon fill="none" stroke="%s" stroke-width="%0.1fpt" stroke-opacity="%0.2f" %s points="',...
                edgecolorname, width, edge_opacity, pattern);
            fprintf(fid,'%0.3f,%0.3f ',[x;y]);
            fprintf(fid,'"/>\n');        
        end
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a line segment
% this algorthm was optimized for large segement counts
function line2svg(fid,group,axpos,x,y,scolorname,style,width)
if ~strcmp(style,'none')
    switch style
        case '--',pattern='stroke-dasharray="8pt,2pt"';
        case ':',pattern='stroke-dasharray="2pt,2pt"';
        case '-.',pattern='stroke-dasharray="8pt,2pt,2pt,2pt"';
        otherwise,pattern='stroke-dasharray="none"';   
    end
    if (isnan(x)==zeros(size(x))&isnan(y)==zeros(size(y)))
        for j=1:20000:length(x)
            xx=x(j:min(length(x),j+19999));
            yy=y(j:min(length(y),j+19999));
            fprintf(fid,'      <polyline fill="none" stroke="%s" stroke-width="%0.1fpt" %s points="',scolorname, width, pattern);
            fprintf(fid,'%0.3f,%0.3f ',[xx;yy]);
            fprintf(fid,'"/>\n');
        end
    else
        parts=find(isnan(x)+isnan(y));
        if parts(1)~=1
            parts=[0 parts];
        end
        if parts(length(parts))~=length(x)
            parts=[parts length(x)+1];
        end
        for j=1:(length(parts)-1)
            xx=x((parts(j)+1):(parts(j+1)-1));
            yy=y((parts(j)+1):(parts(j+1)-1));
            if length(xx)~=0
                fprintf(fid,'      <polyline fill="none" stroke="%s" stroke-width="%0.1fpt" %s points="',scolorname, width, pattern);
                fprintf(fid,'%0.3f,%0.3f ',[xx;yy]);
                fprintf(fid,'"/>\n');
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a circle
function circle2svg(fid,group,axpos,x,y,radius,markeredgecolorname,markerfacecolorname,width)
for j=1:length(x)
    if ~(isnan(x(j)) | isnan(y(j)))
        if ~strcmp(markeredgecolorname,'none') || ~strcmp(markerfacecolorname,'none')
            fprintf(fid,'<circle cx="%0.3f" cy="%0.3f" r="%0.3f" fill="%s" stroke="%s" stroke-width="%0.1fpt" />\n',x(j),y(j),radius,markerfacecolorname,markeredgecolorname,width);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function control2svg(fid,id,ax,group,paperpos)
global PLOT2SVG_globals
set(ax,'Units','pixels');
pos=get(ax,'Position');
pict=getframe(id,pos);
if isempty(pict.colormap)
    pict.colormap=colormap;
end
filename = [PLOT2SVG_globals.basefilename sprintf('%03d',PLOT2SVG_globals.figurenumber) '.' PLOT2SVG_globals.pixelfiletype];
PLOT2SVG_globals.figurenumber = PLOT2SVG_globals.figurenumber + 1;
if isempty(PLOT2SVG_globals.basefilepath)
    current_path = pwd;
else
    current_path = PLOT2SVG_globals.basefilepath;
end
if exist(fullfile(current_path,filename),'file')
    lastwarn('');
    delete(filename);
    if strcmp(lastwarn,'File not found or permission denied.')
        error('Cannot write image file. Make sure that no image is opened in an other program.')    
    end
end
imwrite(pict.cdata,fullfile(PLOT2SVG_globals.basefilepath,filename),PLOT2SVG_globals.pixelfiletype);
set(ax,'Units','normalized');
posNorm=get(ax,'Position');
posInches(1)=posNorm(1)*paperpos(3);
posInches(2)=posNorm(2)*paperpos(4);
posInches(3)=posNorm(3)*paperpos(3);
posInches(4)=posNorm(4)*paperpos(4);
lx = posInches(3);
ly = posInches(4);
pointsx = posInches(1);
pointsy = paperpos(4)-posInches(2)-posInches(4);
fprintf(fid,'<image x="%0.3f" y="%0.3f" width="%0.3f" height="%0.3f" image-rendering="optimizeSpeed" xlink:href="%s" />\n', pointsx, pointsy, lx, ly, filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a text in the axis frame
% the position of the text has to be adapted to the axis scaling
function text2svg(fid,group,axpos,paperpos,id,ax,projection)
originalTextUnits=get(id,'Units');
set(id,'Units','Data');
textpos=get(id,'Position');
set(id,'FontUnits','points');
textfontsize=get(id,'FontSize');
fontsize=convertunit(get(id,'FontSize'),get(id,'FontUnits'),'points');   % convert fontsize to inches
paperposOriginal=get(gcf,'Position');
fontsize=fontsize*paperpos(4)/paperposOriginal(4);
font_color=searchcolor(id,get(id,'Color'));
if strcmp(get(ax,'XScale'),'log')
    textpos(1)=log10(textpos(1));
end
if strcmp(get(ax,'YScale'),'log')
    textpos(2)=log10(textpos(2));
end
if strcmp(get(ax,'ZScale'),'log')
    textpos(3)=log10(textpos(3));
end
[x,y,z] = project(textpos(1),textpos(2),textpos(3),projection);
x = (x*axpos(3)+axpos(1))*paperpos(3);
y = (1-(y*axpos(4)+axpos(2)))*paperpos(4);
textvalign=get(id,'VerticalAlignment');
textalign=get(id,'HorizontalAlignment');
texttext=get(id,'String');
textrot=get(id,'Rotation');
lines=max(size(get(id,'String'),1),1);
if size(texttext,2)~=0
    j=1;
    for i=0:1:(lines-1)
        if iscell(texttext)
            label2svg(fid,group,axpos,id,x,y+i*(fontsize*1.11),convertString(texttext{j}),textalign,textrot,textvalign,lines,paperpos,font_color,0)
        else
            label2svg(fid,group,axpos,id,x,y+i*(fontsize*1.11),convertString(texttext(j,:)),textalign,textrot,textvalign,lines,paperpos,font_color,0)
        end
        j=j+1;   
    end
else
    label2svg(fid,group,axpos,id,x,y,'',textalign,textrot,textvalign,lines,paperpos,font_color,0)
end
set(id,'Units',originalTextUnits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adds the exponents to the axis thickmarks if needed
% MATLAB itself offers no information about this exponent scaling
% the exponent have therefore to be extracted from the thickmarks
function exponent2svg(fid,group,axpos,paperpos,ax)
global PLOT2SVG_globals
if strcmp(get(ax,'XTickLabelMode'),'auto') && strcmp(get(ax,'XScale'),'linear')
    fontsize=convertunit(get(ax,'FontSize'),get(ax,'FontUnits'),'points');   % convert fontsize to inches
    font_color=searchcolor(ax,get(ax,'XColor'));
    if PLOT2SVG_globals.octave
        % Octave stores XTickLabel in a cell array, which does not work nicely with str2num. --Jakob Malm
        axlabelx = get(ax, 'XTickLabel');
        numlabels = zeros(length(axlabelx), 1);
        for ix = 1:length (axlabelx)
            numlabels(ix) = str2num(axlabelx{ix});
        end
    else
        numlabels = str2num(get(ax,'XTickLabel'));
    end
    labelpos=get(ax,'XTick');
    numlabels=numlabels(:);
    labelpos=labelpos(:);
    indexnz=find(labelpos ~= 0);
    if (~isempty(indexnz) && ~isempty(numlabels))
        ratio=numlabels(indexnz)./labelpos(indexnz);
        if round(log10(ratio(1))) ~= 0
            exptext=sprintf('&#215; 10<tspan font-size="%0.1fpt" dy="%0.1fpt">%g</tspan>',0.6*fontsize,-0.6*fontsize,-log10(ratio(1)));
            label2svg(fid,group,axpos,ax,(axpos(1)+axpos(3))*paperpos(3),(1-axpos(2))*paperpos(4)+3*fontsize,exptext,'right',0,'top',1,paperpos,font_color,0)           
        end
    end
end
if strcmp(get(ax,'YTickLabelMode'),'auto') && strcmp(get(ax,'YScale'),'linear')
    fontsize=convertunit(get(ax,'FontSize'),get(ax,'FontUnits'),'points');
    font_color=searchcolor(ax,get(ax,'YColor'));
    if PLOT2SVG_globals.octave
        % Octave stores YTickLabel in a cell array, which does not work nicely with str2num. --Jakob Malm
        axlabely = get(ax, 'YTickLabel');
        numlabels = zeros(length(axlabely), 1);
        for ix = 1:length(axlabely)
            numlabels(ix) = str2num(axlabely{ix});
        end        
    else
        numlabels = str2num(get(ax,'YTickLabel'));
    end
    labelpos=get(ax,'YTick');
    numlabels=numlabels(:);
    labelpos=labelpos(:);
    indexnz=find(labelpos ~= 0);
    if (~isempty(indexnz) && ~isempty(numlabels))
        ratio=numlabels(indexnz)./labelpos(indexnz);
        if round(log10(ratio(1))) ~= 0
            exptext=sprintf('&#215; 10<tspan font-size="%0.1fpt" dy="%0.1fpt">%g</tspan>',0.6*fontsize,-0.6*fontsize,-log10(ratio(1)));
            label2svg(fid,group,axpos,ax,axpos(1)*paperpos(3),(1-(axpos(2)+axpos(4)))*paperpos(4)-0.5*fontsize,exptext,'left',0,'bottom',1,paperpos,font_color,0)           
        end
    end
end
if strcmp(get(ax,'ZTickLabelMode'),'auto') && strcmp(get(ax,'ZScale'),'linear')
    fontsize=convertunit(get(ax,'FontSize'),get(ax,'FontUnits'),'points');
    font_color=searchcolor(ax,get(ax,'ZColor'));
    if PLOT2SVG_globals.octave
        % Octave stores ZTickLabel in a cell array, which does not work nicely with str2num. --Jakob Malm
        axlabelz = get (ax, 'ZTickLabel');
        numlabels = zeros(length(axlabelz), 1);
        for ix = 1:length(axlabelz)
            numlabels(ix) = str2num(axlabelz{ix});
        end
    else
        numlabels = str2num(get(ax,'ZTickLabel'));
    end
    labelpos=get(ax,'ZTick');
    numlabels=numlabels(:);
    labelpos=labelpos(:);
    indexnz=find(labelpos ~= 0);
    if (~isempty(indexnz) && ~isempty(numlabels))
        ratio=numlabels(indexnz)./labelpos(indexnz);
        if round(log10(ratio(1))) ~= 0
            exptext=sprintf('&#215; 10<tspan font-size="%0.1fpt" dy="%0.1fpt">%g</tspan>',0.6*fontsize,-0.6*fontsize,-log10(ratio(1)));
            label2svg(fid,group,axpos,ax,axpos(1)*paperpos(3),(1-(axpos(2)+axpos(4)))*paperpos(4)-0.5*fontsize,exptext,'left',0,'top',1,paperpos,font_color,0)           
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a label in the figure
% former versions of FrameMaker supported the commands FDY and FDX to shift the text
% this commands were replaced by a shift parameter that is normed by the font size
function label2svg(fid,group,axpos,id,x,y,tex,align,angle,valign,lines,paperpos,font_color,exponent)
if isempty(tex)
    return;
end
textfontname=get(id,'FontName');
set(id,'FontUnits','points');
textfontsize=get(id,'FontSize');
if isfield(get(id),'Interpreter')
    if strcmp(get(id,'Interpreter'),'tex')
        latex=1;
    else
        latex=0;
    end
else
    latex=1;
end
fontsize=convertunit(get(id,'FontSize'),get(id,'FontUnits'),'points');   % convert fontsize to inches
paperposOriginal=get(gcf,'Position');
fontsize=fontsize*paperpos(4)/paperposOriginal(4);
textfontsize=textfontsize*paperpos(4)/paperposOriginal(4);
switch lower(valign)
    case 'top',shift=fontsize*0.8;
    case 'cap',shift=fontsize*0.7;
    case 'middle',shift=-((lines-1)/2*fontsize*1.11)+fontsize*0.3;
    case 'bottom',shift=-((lines-1)*fontsize*1.11)+fontsize*-0.04;
    otherwise,shift=0;
end
switch lower(align)
    case 'right', anchor = 'end'; 
    case 'center',anchor = 'middle';
    otherwise,anchor = 'start';
end
if iscellstr(tex)
    tex = strvcat(tex);
elseif ~ ischar(tex)
    error('Invalid character type');
end    
if latex==1 
    tex=strrep(tex,'\alpha','&#945;');
    tex=strrep(tex,'\beta','&#946;');
    tex=strrep(tex,'\gamma','&#947;');
    tex=strrep(tex,'\delta','&#948;');
    tex=strrep(tex,'\epsilon','&#949;');
    tex=strrep(tex,'\zeta','&#950;');
    tex=strrep(tex,'\eta','&#951;');
    tex=strrep(tex,'\theta','&#952;');
    tex=strrep(tex,'\vartheta','&#977;');
    tex=strrep(tex,'\iota','&#953;');
    tex=strrep(tex,'\kappa','&#954;');
    tex=strrep(tex,'\lambda','&#955;');
    tex=strrep(tex,'\mu','&#181;');
    tex=strrep(tex,'\nu','&#957;');
    tex=strrep(tex,'\xi','&#958;');
    tex=strrep(tex,'\pi','&#960;');
    tex=strrep(tex,'\roh','&#961;');
    tex=strrep(tex,'\sigma','&#963;');
    tex=strrep(tex,'\varsigma','&#962;');
    tex=strrep(tex,'\tau','&#964;');
    tex=strrep(tex,'\upsilon','&#965;');
    tex=strrep(tex,'\phi','&#966;');
    tex=strrep(tex,'\chi','&#967;');
    tex=strrep(tex,'\psi','&#968;');
    tex=strrep(tex,'\omega','&#969;');
    tex=strrep(tex,'\Gamma','&#915;');
    tex=strrep(tex,'\Delta','&#916;');
    tex=strrep(tex,'\Theta','&#920;');
    tex=strrep(tex,'\Lambda','&#923;');
    tex=strrep(tex,'\Xi','&#926;');
    tex=strrep(tex,'\Pi','&#928;');
    tex=strrep(tex,'\Sigma','&#931;');
    tex=strrep(tex,'\Tau','&#932;');
    tex=strrep(tex,'\Upsilon','&#933;');
    tex=strrep(tex,'\Phi','&#934;');
    tex=strrep(tex,'\Psi','&#936;');
    tex=strrep(tex,'\Omega','&#937;');
    tex=strrep(tex,'\infty','&#8734;');
    tex=strrep(tex,'\pm','&#177;');
    tex=strrep(tex,'\Im','&#8465;');
    tex=strrep(tex,'\Re','&#8476;');
    tex=strrep(tex,'\approx','&#8773;');
    tex=strrep(tex,'\leq','&#8804;');
    tex=strrep(tex,'\geq','&#8805;');
    tex=strrep(tex,'\times','&#215;');
    tex=strrep(tex,'\leftrightarrow','&#8596;');
    tex=strrep(tex,'\leftarrow','&#8592;');
    tex=strrep(tex,'\uparrow','&#8593;');
    tex=strrep(tex,'\rightarrow','&#8594;');
    tex=strrep(tex,'\downarrow','&#8595;');
    tex=strrep(tex,'\circ','&#186;');
    tex=strrep(tex,'\propto','&#8733;');
    tex=strrep(tex,'\partial','&#8706;');
    tex=strrep(tex,'\bullet','&#8226;');
    tex=strrep(tex,'\div','&#247;');
    tex=latex2svg(tex,textfontname,textfontsize,'FNormal');
end
if isempty(tex)
    return;
end
if exponent
    tex=sprintf('10<tspan font-size="%0.1fpt" dy="%0.1fpt">%s</tspan>',0.6*textfontsize,-0.6*textfontsize,tex);
    shift = shift + 0.4*fontsize;   % Small correction to make it look nicer
end
fprintf(fid,'  <g transform="translate(%0.3f,%0.3f)">\n',x,y+shift);
fprintf(fid,'    <g transform="rotate(%0.1f)">\n',-angle);
fprintf(fid,'      <text x="%0.3f" y="%0.3f" font-family="%s" text-anchor="%s" font-size="%0.1fpt" fill="%s" >', 0, 0, textfontname, anchor, textfontsize * 0.8, font_color);
fprintf(fid,'%s',tex);
fprintf(fid,'</text>\n'); 
fprintf(fid,'    </g>\n');
fprintf(fid,'  </g>\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% converts LATEX strings into Framemaker strings
function returnvalue=latex2svg(StringText,font,size,style)
if isempty(StringText)
    returnvalue='';
else
    leftbracket=0;
    rightbracket=0;
    bracketcounter=0;
    leftbracketpos=[];
    rightbracketpos=[];
    returnvalue=[];
    for i=1:length(StringText)
        if rightbracket==leftbracket
            returnvalue=[returnvalue StringText(i)];    
        end
        if StringText(i)=='{'
            leftbracket=leftbracket+1;
            bracketcounter=bracketcounter+1;
            leftbracketpos=[leftbracketpos i];
        end
        if StringText(i)=='}'
            rightbracket=rightbracket+1;
            rightbracketpos=[rightbracketpos i];
            if rightbracket==leftbracket
                fontnew=font;
                sizenew=size;
                stylenew=style;
                if leftbracketpos(leftbracket-bracketcounter+1)~=1
                    switch StringText(leftbracketpos(leftbracket-bracketcounter+1)-1)   
                        case '^'
                            stylenew='super';
                            returnvalue=returnvalue(1:(end-1));
                        case '_'
                            stylenew='sub';
                            returnvalue=returnvalue(1:(end-1));
                    end
                end
                if strcmp(style,stylenew)
                    format=[];
                    formatend=[];
                else
                    format=['<tspan baseline-shift="' stylenew '">'];
                    formatend='</tspan>';
                end
                textinbrackets=StringText((leftbracketpos(leftbracket-bracketcounter+1)+1):(rightbracketpos(rightbracket)-1));
                foundpos=findstr(textinbrackets,'\bf');
                if ~isempty(foundpos)
                    textinbrackets=strrep(textinbrackets,'\bf','<tspan font-weight="bold">');
                    textinbrackets=[textinbrackets '</tspan>'];
                end
                foundpos=findstr(textinbrackets,'\it');
                if ~isempty(foundpos)
                    textinbrackets=strrep(textinbrackets,'\it','<tspan font-style="italic">');
                    textinbrackets=[textinbrackets '</tspan>'];
                end
                returnvalue=[returnvalue(1:(end-1)) format latex2svg(textinbrackets,fontnew,sizenew,stylenew) formatend];
                bracketcounter=0;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function name=searchcolor(id,value)
if ischar(value)
    name = value;
else
    name=sprintf('#%02x%02x%02x',fix(value(1)*255),fix(value(2)*255),fix(value(3)*255));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rvalue=convertunit(value,from,to)
switch lower(from)  % convert from input unit to points
    case 'points', rvalue=value;
    case 'centimeters', rvalue=value/2.54*72;
    case 'inches', rvalue=value*72; % 72 points = 1 inch
    otherwise, error(['Unknown unit ' from '.']);
end
switch lower(to)    % convert from points to specified unit
    case 'points'; % do nothing
    case 'centimeters', rvalue=rvalue*2.54/72;
    case 'inches', rvalue=rvalue/72;    % 72 points = 1 inch
    otherwise, error(['Unknown unit ' to '.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function strString=addBackSlash( strSlash)
% adds a backslash at the last position of the string (if not already there)
if ( strSlash(end) ~= filesep)
    strString = [ strSlash filesep];
else
    strString = strSlash;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function strExt=getFileExtension( strFileName)
% returns the file extension of a filename
[path, name, strExt] = fileparts( strFileName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function StringText=convertString(StringText)
if iscell(StringText) % Octave stores some strings in cell arrays. --Jakob Malm
    StringText = StringText{1};
end
if ~isempty(StringText)
    StringText=strrep(StringText,'&','&amp;');  % Do not change sequence !!
    StringText=strrep(StringText,'\\','\');
    StringText=strrep(StringText,'<','&lt;');
    StringText=strrep(StringText,'>','&gt;');
    StringText=strrep(StringText,'"','&quot;');
    % Workaround for Firefox and Inkscape
    StringText=strrep(StringText,'°','Â°');
    %StringText=strrep(StringText,'°','&deg;');
    StringText=strrep(StringText,'±','&plusmn;');
    StringText=strrep(StringText,'µ','&micro;');
    StringText=strrep(StringText,'²','&sup2;');
    StringText=strrep(StringText,'³','&sup3;');
    StringText=strrep(StringText,'¼','&frac14;');
    StringText=strrep(StringText,'½','&frac12;');
    StringText=strrep(StringText,'¾','&frac34;');
    StringText=strrep(StringText,'©','&copy;');
    StringText=strrep(StringText,'®','&reg;');
    StringText=deblank(StringText);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IdString = createId
global PLOT2SVG_globals
IdString = ['ID' sprintf('%06d',PLOT2SVG_globals.runningIdNumber)];
PLOT2SVG_globals.runningIdNumber = PLOT2SVG_globals.runningIdNumber + 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [projection, edges] = get_projection(ax,id)
global PLOT2SVG_globals
xc = get(ax,'CameraTarget');
phi = get(ax,'CameraViewAngle');
vi = get(ax,'View');
xi = get(ax,'XLim');
yi = get(ax,'YLim');
zi = get(ax,'ZLim');
if strcmp(get(ax,'XScale'),'log')
    xi=log10(get(ax,'XLim'));
end
if strcmp(get(ax,'YScale'),'log')
    yi=log10(get(ax,'YLim'));
end
if strcmp(get(ax,'ZScale'),'log')
    zi=log10(get(ax,'ZLim'));
end
projection.xi = xi;
projection.yi = yi;
projection.zi = zi;
projection.aspect_scaling = get(ax,'DataAspectRatio');
xc(1) = (xc(1) - xi(1))/(xi(2)-xi(1));
xc(2) = (xc(2) - yi(1))/(yi(2)-yi(1));
xc(3) = (xc(3) - zi(1))/(zi(2)-zi(1));
x = [xi(1) xi(2) xi(1) xi(2) xi(1) xi(2) xi(1) xi(2)]/projection.aspect_scaling(1);
y = [yi(1) yi(1) yi(2) yi(2) yi(1) yi(1) yi(2) yi(2)]/projection.aspect_scaling(2);
z = [zi(1) zi(1) zi(1) zi(1) zi(2) zi(2) zi(2) zi(2)]/projection.aspect_scaling(3);
if PLOT2SVG_globals.octave
        projection.A = get(ax,'x_ViewTransform');
        projection.A(3,:) = -projection.A(3,:);
        projection.A(1:3,4) = 0;
else
    if strcmp(get(ax,'Projection'),'orthographic')
        projection.A = viewmtx(vi(1),vi(2));
    else
        projection.A = viewmtx(vi(1),vi(2),phi,xc);
    end
end
if (vi(1) == 0) && (mod(vi(2),90) == 0)
    projection.xyplane = true;
else
    projection.xyplane = false;
end
axpos = get(ax,'Position');
figpos = get(id,'Position');
[m,n] = size(x);
x4d = [x(:),y(:),z(:),ones(m*n,1)]';
x2d = projection.A*x4d;
x2 = zeros(m,n); y2 = zeros(m,n); z2 = zeros(m,n);
x2(:) = x2d(1,:)./x2d(4,:);
y2(:) = x2d(2,:)./x2d(4,:);
projection.ax = ax;
projection.xrange = max(x2) - min(x2);
projection.yrange = max(y2) - min(y2);
projection.xoffset = (max(x2) + min(x2))/2;
projection.yoffset = (max(y2) + min(y2))/2;
if (strcmp(get(ax,'PlotBoxAspectRatioMode'),'manual') || strcmp(get(ax,'DataAspectRatioMode'),'manual'))
    if (projection.xrange*axpos(4)*figpos(4) < projection.yrange*axpos(3)*figpos(3))
        projection.xrange = projection.yrange*axpos(3)*figpos(3)/axpos(4)/figpos(4);
    else
        projection.yrange = projection.xrange*axpos(4)*figpos(4)/axpos(3)/figpos(3);
    end
end
x2(:) = (x2d(1,:)./x2d(4,:) - projection.xoffset)/projection.xrange + 0.5;
y2(:) = (x2d(2,:)./x2d(4,:) - projection.yoffset)/projection.yrange + 0.5;
z2(:) =  x2d(3,:);
edges = [x2; y2; z2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x2,y2,z2] = project(x,y,z,projection)
[m,n] = size(x);
if strcmp(get(projection.ax,'XDir'),'reverse')
    xi = projection.xi;
    x = (1 - (x - xi(1)) / (xi(2) - xi(1))) * (xi(2) - xi(1)) + xi(1);
end
if strcmp(get(projection.ax,'YDir'),'reverse')
    yi = projection.yi;
    y = (1 - (y - yi(1)) / (yi(2) - yi(1))) * (yi(2) - yi(1)) + yi(1);
end
x4d = [x(:)/projection.aspect_scaling(1), y(:)/projection.aspect_scaling(2), z(:)/projection.aspect_scaling(3), ones(m*n,1)]';
x2d = projection.A*x4d;
x2 = zeros(m,n); y2 = zeros(m,n); z2 = zeros(m,n);
x2(:) = (x2d(1,:)./x2d(4,:) - projection.xoffset)/projection.xrange + 0.5;
y2(:) = (x2d(2,:)./x2d(4,:) - projection.yoffset)/projection.yrange + 0.5;
z2(:) =  x2d(3,:);
%x = [0 1 0 1 0 1 0 1];
%y = [0 0 1 1 0 0 1 1];
%z = [0 0 0 0 1 1 1 1];


function [f, v, fvc, fva] = surface2patch(s)
x = get(s, 'xdata');
y = get(s, 'ydata');
z = get(s, 'zdata');
c = get(s, 'cdata');
a = get(s, 'AlphaData');
if ~isempty(x) & ~isequal(size(x),size(z))
    x = repmat(x(:)',size(z,1),1);
end
if ~isempty(y) & ~isequal(size(y),size(z))
  y = repmat(y(:),1,size(z,2));
end
[m n]= size(z);
if isempty(x)
    [x y] = meshgrid(1:n, 1:m);
end
[cm cn cp] = size(c);
[am an ap] = size(a);
%if cm==(m-1) & cn==(n-1)
%    cmode = 'f';
%elseif cm==m & cn==n
%    cmode = 'v';
%else  
%    cmode = '';
%end
v = [x(:) y(:) z(:)];
q = [1:m*n-m-1]';
q(m:m:end) = [];
fvc = reshape(c, [cm*cn cp]);
fva = reshape(a, [am*an ap]);
f = [q q+m q+m+1 q+1];
