function [ Label bRefresh bCancel ] = ExpertUserInterface(cCentroid, cCloserExamples, cDistantExamples)

Ventricular_ptr = @(a,b)(ChooseLabel(a,b,'Ventricular'));
Normal_ptr = @(a,b)(ChooseLabel(a,b,'Normal'));
Supraventricular_ptr = @(a,b)(ChooseLabel(a,b,'Supraventricular'));
Unknown_ptr = @(a,b)(ChooseLabel(a,b,'Unknown'));

% open a new figure
fig_hnd = figure();

set(fig_hnd, 'Tag', 'a2hbc' );
screen_size = get(0,'ScreenSize');
set(fig_hnd, 'Name', 'Expert interface' );
set(fig_hnd, 'position', [10 screen_size(4)*0.43 screen_size(3)*0.97 screen_size(4)*0.5] );

%no expert opinion
user_data.Label = [];
user_data.Refresh = false;
user_data.Cancel = false;

%% Plots
CentroidPanel = uipanel(    'Title', [ 'Centroid heartbeat. Total of ' num2str(cCentroid{7}) ' heartbeats in this cluster.' ], ... 
                            'FontSize',12,...
                            'Position', [0.0011    0.0177    0.3146    0.9693]);

ECG_centroid_hdl = axes('Parent', CentroidPanel);
set(ECG_centroid_hdl, 'units','normalized' );
set(ECG_centroid_hdl, 'position', [0.08    0.48    0.9096    0.4638] )

RR_centroid_hdl = axes('Parent', CentroidPanel);
set(RR_centroid_hdl, 'units','normalized' );
set(RR_centroid_hdl, 'position', [0.08    0.045    0.36    0.35] )

Morph_centroid_hdl = axes('Parent', CentroidPanel);
set(Morph_centroid_hdl, 'units','normalized' );
set(Morph_centroid_hdl, 'position', [0.44    0.045    0.5402    0.35] )


CloserPanel = uipanel(    'Title', [ 'Closer heartbeats. Showing  ' num2str(length(cCloserExamples{1})) ' heartbeats.' ], ... 
                            'FontSize',12,...
                            'Position', [0.3171    0.0177    0.3146    0.9693]);

ECG_closer_hdl = axes('Parent', CloserPanel);
set(ECG_closer_hdl, 'units','normalized' );
set(ECG_closer_hdl, 'position', [0.08    0.48    0.9096    0.4638] )

RR_closer_hdl = axes('Parent', CloserPanel);
set(RR_closer_hdl, 'units','normalized' );
set(RR_closer_hdl, 'position', [0.08    0.045    0.36    0.35] )

Morph_closer_hdl = axes('Parent', CloserPanel);
set(Morph_closer_hdl, 'units','normalized' );
set(Morph_closer_hdl, 'position', [0.44    0.045    0.5402    0.35] )

DistantPanel = uipanel(    'Title',[ 'Distant heartbeats. Showing  ' num2str(length(cDistantExamples{1})) ' heartbeats.' ], ... 
                            'FontSize',12,...
                            'Position', [0.6336    0.0177    0.3146    0.9693]);

ECG_distant_hdl = axes('Parent', DistantPanel);
set(ECG_distant_hdl, 'units','normalized' );
set(ECG_distant_hdl, 'position', [0.08    0.48    0.9096    0.4638] )

RR_distant_hdl = axes('Parent', DistantPanel);
set(RR_distant_hdl, 'units','normalized' );
set(RR_distant_hdl, 'position', [0.08    0.045    0.36    0.35] )

Morph_distant_hdl = axes('Parent', DistantPanel);
set(Morph_distant_hdl, 'units','normalized' );
set(Morph_distant_hdl, 'position', [0.44    0.045    0.5402    0.35] )

%% Buttons

ControlPanel = uipanel( 'Position', [0.9503    0.0177    0.0485    0.9480]);
                      
uicontrol('style','pushbutton', ...
          'Parent', ControlPanel, ...
          'string', 'Ventricular', ...
          'FontName','Arial', ...
          'FontSize',10, ...
          'FontWeight','Bold', ...
          'ForegroundColor',[1 0 0], ...
          'units','normalized', ...
          'position', [ 0.05    0.7119    0.9    0.0500 ], ...
          'callback', Ventricular_ptr );

uicontrol('style','pushbutton', ...
          'Parent', ControlPanel, ...
          'string', 'Supraventricular', ...
          'FontName','Arial', ...
          'FontSize',10, ...
          'FontWeight','Bold', ...
          'ForegroundColor',[0 0.5 0], ...
          'units','normalized', ...
          'position', [ 0.05    0.7843    0.9    0.0500 ], ...
          'callback', Supraventricular_ptr);

uicontrol('style','pushbutton', ...
          'Parent', ControlPanel, ...
          'string', 'Normal', ...
          'FontName','Arial', ...
          'FontSize',10, ...
          'FontWeight','Bold', ...
          'ForegroundColor',[0 0 1], ...
          'units','normalized', ...
          'position', [ 0.05    0.8586    0.9    0.0500 ], ...
          'callback', Normal_ptr);

uicontrol('style','pushbutton', ...
          'Parent', ControlPanel, ...
          'string', 'Unknown', ...
          'FontName','Arial', ...
          'FontSize',10, ...
          'FontWeight','Bold', ...
          'ForegroundColor',[1 1 1], ...
          'units','normalized', ...
          'position', [ 0.05    0.64    0.9    0.0500 ], ...
          'callback', Unknown_ptr);

uicontrol('style','pushbutton', ...
          'Parent', ControlPanel, ...
          'string', 'Skip', ...
          'FontName','Arial', ...
          'FontSize',10, ...
          'FontWeight','Bold', ...
          'ForegroundColor',[0 0 1], ...
          'units','normalized', ...
          'position', [ 0.05    0.5    0.9    0.0500 ], ...
          'callback', @SkipCluster);

uicontrol('style','pushbutton', ...
          'Parent', ControlPanel, ...
          'string', 'Refresh', ...
          'FontName','Arial', ...
          'FontSize',10, ...
          'FontWeight','Bold', ...
          'ForegroundColor',[0 0 1], ...
          'units','normalized', ...
          'position', [ 0.05    0.43    0.9    0.0500 ], ...
          'callback',@RefreshExamples);

uicontrol('style','pushbutton', ...
          'Parent', ControlPanel, ...
          'string', 'Cancel', ...
          'FontName','Arial', ...
          'FontSize',10, ...
          'FontWeight','Bold', ...
          'ForegroundColor',[1 0 0], ...
          'units','normalized', ...
          'position', [ 0.05    0.05    0.9    0.0500 ], ...
          'callback',@CancelButton);

      
%% plot data

set(fig_hnd,'ToolBar','figure');
    
%plot signals

% cCloserExamples = { aux_rel_posQRS ... %relative position of the heartbeat
%                     aux_ECGslices ... % ecg excerpt
%                     aux_rel_posRR ... % relative position of the heartbeat 
%                     aux_RRslices ... % RR interval sequence
%                     };
% cCentroid, cCloserExamples, cDistantExamples

box(ECG_closer_hdl, 'off');
box(RR_closer_hdl, 'off');
box(ECG_distant_hdl, 'off');
box(RR_distant_hdl, 'off');

%% Centroid heartbeat

plot( ECG_centroid_hdl,  cCentroid{2} );
hold(ECG_centroid_hdl, 'on');
plot( ECG_centroid_hdl,  [ cCentroid{1} cCentroid{1}], [ min(min(cCentroid{2})) max(max(cCentroid{2}))], '--r' );
hold(ECG_centroid_hdl, 'off');
TightAxes(ECG_centroid_hdl);
TimeRes = 0.4;
TimeGrid = fliplr(cCentroid{1}:-TimeRes*360:1);
TimeGrid = [ TimeGrid(1:end-1) cCentroid{1}:TimeRes*360:size(cCentroid{2},1)];
set(ECG_centroid_hdl, 'XTick', colvec(TimeGrid));
set(ECG_centroid_hdl, 'XTickLabel', Seconds2HMS((TimeGrid-cCentroid{1})/360) );
title(ECG_centroid_hdl, 'First Two PCA Components' );
box(ECG_centroid_hdl, 'off');

plot( RR_centroid_hdl,  cCentroid{4}, '--ob' );
hold(RR_centroid_hdl, 'on');
plot( RR_centroid_hdl,  [ cCentroid{3} cCentroid{3}], [ min(min(cCentroid{4})) max(max(cCentroid{4}))], '--r' );
hold(RR_centroid_hdl, 'off');
TightAxes(RR_centroid_hdl);
set(RR_centroid_hdl, 'XTick', cCentroid{3});
set(RR_centroid_hdl, 'XTickLabel', 'Current heartbeat' );
title(RR_centroid_hdl, 'Local RR interval evolution' );
box(RR_centroid_hdl, 'off');

plot( Morph_centroid_hdl,  cCentroid{6} );
hold(Morph_centroid_hdl, 'on');
plot( Morph_centroid_hdl,  [ cCentroid{5} cCentroid{5}], [ min(min(cCentroid{6})) max(max(cCentroid{6}))], '--r' );
hold(Morph_centroid_hdl, 'off');
TightAxes(Morph_centroid_hdl);
TimeRes = 0.1;
TimeGrid = fliplr(cCentroid{5}:-TimeRes*360:1);
TimeGrid = [ TimeGrid(1:end-1) cCentroid{5}:TimeRes*360:size(cCentroid{6},1)];
set(Morph_centroid_hdl, 'XTick', colvec(TimeGrid));
set(Morph_centroid_hdl, 'YTick', []);
set(Morph_centroid_hdl, 'XTickLabel', Seconds2HMS((TimeGrid-cCentroid{5})/360) );
title(Morph_centroid_hdl, 'Morphology details' );
box(Morph_centroid_hdl, 'off');


%% Closer heartbeats

QRS_positions = cCloserExamples{1};
ECG = cCloserExamples{2};

plot( ECG_closer_hdl, (1:size(ECG{1},1))-QRS_positions(1), ECG{1} );
hold(ECG_closer_hdl, 'on');
cellfun(@(a,b)(plot( ECG_closer_hdl, (1:size(a,1))-b, a )), ECG(min(end,2):end), num2cell(QRS_positions(min(end,2):end)) );
plot( ECG_closer_hdl,  [ 0 0 ], [ min(min(cell2mat(ECG))) max(max(cell2mat(ECG)))], '--r' );
hold(ECG_closer_hdl, 'off');
TightAxes(ECG_closer_hdl);
TimeRes = 0.4;
TimeGrid = fliplr(QRS_positions(1):-TimeRes*360:1);
TimeGrid = [ TimeGrid(1:end-1) QRS_positions(1):TimeRes*360:size(ECG{1},1)]-QRS_positions(1);
set(ECG_closer_hdl, 'XTick', colvec(TimeGrid));
set(ECG_closer_hdl, 'XTickLabel', Seconds2HMS((TimeGrid)/360) );
title(ECG_closer_hdl, 'First Two PCA Components' );
box(ECG_closer_hdl, 'off');

QRS_positions = cCloserExamples{3};
RR_series = cCloserExamples{4};

plot( RR_closer_hdl, (1:size(RR_series{1},1))-QRS_positions(1), RR_series{1}, ':xb' );
hold(RR_closer_hdl, 'on');
cellfun(@(a,b)(plot( RR_closer_hdl, (1:size(a,1))-b, a, ':xb' )), RR_series(min(end,2):end), num2cell(QRS_positions(min(end,2):end)) );
plot( RR_closer_hdl,  [ 0 0 ], [ min(min(cell2mat(RR_series))) max(max(cell2mat(RR_series)))], '--r' );
hold(RR_closer_hdl, 'off');
TightAxes(RR_closer_hdl);
set(RR_closer_hdl, 'XTick', 0);
set(RR_closer_hdl, 'XTickLabel', 'Current heartbeat' );
title(RR_closer_hdl, 'Local RR interval evolution' );
box(RR_closer_hdl, 'off');

QRS_positions = cCloserExamples{5};
ECG = cCloserExamples{6};

plot( Morph_closer_hdl, (1:size(ECG{1},1))-QRS_positions(1), ECG{1} );
hold(Morph_closer_hdl, 'on');
cellfun(@(a,b)(plot( Morph_closer_hdl, (1:size(a,1))-b, a )), ECG(min(end,2):end), num2cell(QRS_positions(min(end,2):end)) );
hold(Morph_closer_hdl, 'off');
TightAxes(Morph_closer_hdl);
TimeRes = 0.1;
TimeGrid = fliplr(QRS_positions(1):-TimeRes*360:1);
TimeGrid = [ TimeGrid(1:end-1) QRS_positions(1):TimeRes*360:size(ECG{1},1)]-QRS_positions(1);
set(Morph_closer_hdl, 'XTick', colvec(TimeGrid));
set(Morph_closer_hdl, 'XTickLabel', Seconds2HMS((TimeGrid)/360) );
set(Morph_closer_hdl, 'YTick', []);
title(Morph_closer_hdl, 'Morphology details' );
box(Morph_closer_hdl, 'off');

%% Distant heartbeats

QRS_positions = cDistantExamples{1};
ECG = cDistantExamples{2};

plot( ECG_distant_hdl, (1:size(ECG{1},1))-QRS_positions(1), ECG{1} );
hold(ECG_distant_hdl, 'on');
cellfun(@(a,b)(plot( ECG_distant_hdl, (1:size(a,1))-b, a )), ECG(min(end,2):end), num2cell(QRS_positions(min(end,2):end)) );
plot( ECG_distant_hdl,  [ 0 0 ], [ min(min(cell2mat(ECG))) max(max(cell2mat(ECG)))], '--r' );
hold(ECG_distant_hdl, 'off');
TightAxes(ECG_distant_hdl);
TimeRes = 0.4;
TimeGrid = fliplr(QRS_positions(1):-TimeRes*360:1);
TimeGrid = [ TimeGrid(1:end-1) QRS_positions(1):TimeRes*360:size(ECG{1},1)]-QRS_positions(1);
set(ECG_distant_hdl, 'XTick', colvec(TimeGrid));
set(ECG_distant_hdl, 'XTickLabel', Seconds2HMS((TimeGrid)/360) );
title(ECG_distant_hdl, 'First Two PCA Components' );
box(ECG_distant_hdl, 'off');

QRS_positions = cDistantExamples{3};
RR_series = cDistantExamples{4};

plot( RR_distant_hdl, (1:size(RR_series{1},1))-QRS_positions(1), RR_series{1}, ':xb' );
hold(RR_distant_hdl, 'on');
cellfun(@(a,b)(plot( RR_distant_hdl, (1:size(a,1))-b, a, ':xb' )), RR_series(min(end,2):end), num2cell(QRS_positions(min(end,2):end)) );
plot( RR_distant_hdl,  [ 0 0 ], [ min(min(cell2mat(RR_series))) max(max(cell2mat(RR_series)))], '--r' );
hold(RR_distant_hdl, 'off');
TightAxes(RR_distant_hdl);
set(RR_distant_hdl, 'XTick', 0);
set(RR_distant_hdl, 'XTickLabel', 'Current heartbeat' );
title(RR_distant_hdl, 'Local RR interval evolution' );
box(RR_distant_hdl, 'off');

QRS_positions = cDistantExamples{5};
ECG = cDistantExamples{6};

plot( Morph_distant_hdl, (1:size(ECG{1},1))-QRS_positions(1), ECG{1} );
hold(Morph_distant_hdl, 'on');
cellfun(@(a,b)(plot( Morph_distant_hdl, (1:size(a,1))-b, a )), ECG(min(end,2):end), num2cell(QRS_positions(min(end,2):end)) );
hold(Morph_distant_hdl, 'off');
TightAxes(Morph_distant_hdl);
TimeRes = 0.1;
TimeGrid = fliplr(QRS_positions(1):-TimeRes*360:1);
TimeGrid = [ TimeGrid(1:end-1) QRS_positions(1):TimeRes*360:size(ECG{1},1)]-QRS_positions(1);
set(Morph_distant_hdl, 'XTick', colvec(TimeGrid));
set(Morph_distant_hdl, 'XTickLabel', Seconds2HMS((TimeGrid)/360) );
set(Morph_distant_hdl, 'YTick', []);
title(Morph_distant_hdl, 'Morphology details' );
box(Morph_distant_hdl, 'off');


%% Wait for user interaction

set(fig_hnd , 'User', user_data );

uiwait(fig_hnd); 

if( ishandle(fig_hnd) )
    user_data = get(fig_hnd , 'User');

    Label  = user_data.Label;
    bRefresh = user_data.Refresh;
    bCancel = user_data.Cancel;
    
    close(fig_hnd);
else
    bCancel = true;
    Label  = [];
    bRefresh = [];
end


function ChooseLabel(obj, event_obj, strLabel)

fig_hnd = gcf();
user_data = get(fig_hnd , 'User');
user_data.Label = strLabel;
set(fig_hnd , 'User', user_data );
uiresume(fig_hnd); 

function SkipCluster(obj, event_obj)

fig_hnd = gcf();
uiresume(fig_hnd); 

function CancelButton(obj, event_obj)
fig_hnd = gcf();
user_data = get(fig_hnd , 'User');
user_data.Cancel = true;
set(fig_hnd , 'User', user_data );
uiresume(fig_hnd); 

function RefreshExamples(obj, event_obj)

fig_hnd = gcf();
user_data = get(fig_hnd , 'User');
user_data.Refresh = true;
set(fig_hnd , 'User', user_data );
uiresume(fig_hnd); 


function TightAxes(axes_hdl)

axis(axes_hdl, 'tight')
xlimits = xlim(axes_hdl);
xrange = diff(xlimits);
ylimits = ylim(axes_hdl);
yrange = diff(ylimits);
xlim(axes_hdl,  xlimits + [ -0.05*xrange 0.05*xrange ] )
ylim(axes_hdl, ylimits + [ -0.05*yrange 0.05*yrange ] )
