
if( ~exist('UCP_struct', 'var') || ~ishandle(UCP_struct.fig_hdl) )
    UserControlPanel;
else
%     set(UCP_struct.Axes_hdl, 'Visible','off' );
    DisableControPanel;
end

fprintf(1, 'Now you can change controls freely.\n');

% set(UCP_struct.fig_hdl, 'Visible','on' );
EnableControPanel;
uiwait(UCP_struct.fig_hdl);

if( ~ishandle(UCP_struct.fig_hdl) )
    bUserExit = true;
else
%     set(UCP_struct.fig_hdl, 'Visible','off' );
    DisableControPanel;
    
    ParseUserControlPanelInput;
    
end

