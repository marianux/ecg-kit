function DoHouseKeeping(added_paths, previous_open_figure_handles)

our_fig_hdl = [findall(0, 'Tag', 'a2hbc');findall(0, 'Tag', 'progress_bar')];
close(our_fig_hdl);

on_exit_open_handles = findall(0);
handles2_close = setdiff(on_exit_open_handles, previous_open_figure_handles);
delete(handles2_close);

prwaitbar('off')        

rmpath(added_paths);
