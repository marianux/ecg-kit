%% (Internal) function reads the header of signal files
%
% This function creates a report of the signals handled by the ECG wrapper
% object ECG_w. The report includes several views of the signals, from a
% wide to a narrow scale. Some aspects of the report can be configured as
% the detail degree, the length of each time scale and the report format.
% 
% Arguments:
%     
%     +ECG_w: [ECGwrapper object] REQUIRED. The signal handler.
% 
%     +detailLevel: [char] OPTIONAL. The report detail level:
%                   'HighDetail', 'MediumDetail', 'LowDetail'. A higher
%                   detail level means report the whole recording at every
%                   time resolution defined in "win_lengths". High
%                   resolution also means larger reports. Default 'LowDetail'. 
% 
%     +report_mode: [char] OPTIONAL. Information from other tasks like QRS
%                   detection/delineation/classification added to the
%                   signals in case available mode. Possible values are:
%                   'full' 'ECG only' 'QRS detection' 'Wave delineation'
%                   'Heartbeat classification' 
%                   Default "ECG only". 
% 
%     +win_lengths: [numeric] OPTIONAL. The amount and size (in seconds) of
%                   each scale length present in the report. Default 
%                   [60*60 30*60 60 10] it means 1 hour - 30 min - 1 min
%                   and 10 seconds.
% 
%     +report_format: [char] OPTIONAL. The report format of the document.
%                  Default PDF.
% 
%     +filename: [char] OPTIONAL. The report filename. Default same folder
%                of the original recording handled by ECG_w.
% 
% 
% See also plot_ecg_strip
% 
% Author: Mariano Llamedo Soria (llamedom at frba.utn.edu.ar)
% Version: 0.1 beta
% Last update: 17/6/2014
% Birthdate  : 15/8/2012
% Copyright 2008-2015
% 
function reportECG(ECG_w, detailLevel, report_mode, win_lengths, report_format, filename)


%% constants

cKnownModes = {'full' 'ECG only' 'QRS detection' 'Wave delineation' 'Heartbeat classification' };
cKnownFormats = {'pdf' 'eps' 'png' 'tiff' 'jpg' 'bmp' };
cDetails = {'HighDetail' 'MediumDetail' 'LowDetail'};

%% argument parsing

if( nargin < 5 || ~any(strcmpi(cKnownFormats, report_format)) )
    report_format = 'pdf';
end

heasig = ECG_w.ECG_header;

if( nargin < 3 || isempty(intersect(cKnownModes, report_mode)) )
    report_mode = 'ECG only';
end

file_suffix = rowvec(colvec(char(report_mode)'));
prev_length = length(file_suffix);
file_suffix = strrep(file_suffix, '  ', ' ');
while( length(file_suffix) ~= prev_length  )
    prev_length = length(file_suffix);
    file_suffix = strrep(file_suffix, '  ', ' ');
end
file_suffix = ['_' file_suffix];

if( any(strcmpi(report_mode, 'full')) )
    report_mode = cKnownModes(2:end);
end

if( nargin < 6 || isempty(filename) )

%     [report_path, file_recname ]= fileparts(ECG_w.recording_name);
%     report_path = [report_path filesep];
    [~, file_recname ]= fileparts(ECG_w.recording_name);
    report_path = ECG_w.output_path;
    
    if( isfield(heasig, 'recname') )
        filename = [report_path heasig.recname file_suffix '.' report_format];
    else
        filename = [report_path file_recname file_suffix '.' report_format];
    end

    ii = 1;
    while( exist(filename, 'file') )

        if( isfield(heasig, 'recname') )
            filename = [report_path heasig.recname file_suffix '_copy' num2str(ii) '.' report_format];
        else
            filename = [report_path file_recname file_suffix '_copy' num2str(ii) '.' report_format];
        end
        ii = ii + 1;
    end
    
end

if( nargin < 2 || ~any(strcmpi(cDetails, detailLevel)) )
    detailLevel = 'LowDetail'; % seconds
end

if( nargin < 4 || isempty(win_lengths) || ~isnumeric(win_lengths) )
    win_lengths = [60*60 30*60 59 7]; % seconds
end

if( ~isobject(ECG_w) || ~isa(ECG_w, 'ECGwrapper') )
    error('reportECG:BadArg', 'Need an ECGwrapper object to work.')
end

init_ghostscript();

bAux = (win_lengths * heasig.freq) <= heasig.nsamp;

% Activate the progress_struct bar.
pb = progress_bar('Report ECG');

pb.Loops2Do = sum(bAux);

ii = 1;
fig_hdl = [];


try 
    
    for win_length = win_lengths(bAux)

        pb.start_loop();

        if( (win_length * heasig.freq) > heasig.nsamp )
            win_length = round(heasig.nsamp / heasig.freq);
        end

        switch detailLevel
            case 'LowDetail'
                max_cant_figs = 10;
            case 'MediumDetail'
                max_cant_figs = 30;
            case 'HighDetail'
                max_cant_figs = realmax;
        end

        cant_figs = min(max_cant_figs, ceil(heasig.nsamp / win_length / heasig.freq));
        starts = linspace(0, round(heasig.nsamp / heasig.freq) - win_length, cant_figs);

        str_aux = ['Reporting for window size: ' Seconds2HMS(win_length) ];
        pb1 = pb.AddInerLoop( str_aux );
        pb1.Loops2Do = cant_figs;

        jj = 1;

        for start_time = starts

            pb1.start_loop();

            fig_hdl = prepare_fig_hdl( fig_hdl );

            if( all(strcmpi(report_mode, 'ECG only')) )
                plot_ecg_strip(ECG_w, 'PrettyPrint', true, 'Figure_handle', fig_hdl, 'Start_time', start_time, 'End_time', start_time + win_length);
            else

                plot_ecg_strip(ECG_w, 'PrettyPrint', true, ... 
                                      'Figure_handle', fig_hdl, ...
                                      'Start_time', start_time, ...
                                      'End_time', start_time + win_length);
            end

            set(fig_hdl, 'Visible', 'off');
            str_aux2 = [ str_aux ' - Start: ' Seconds2HMS(start_time) ];
            pb1.checkpoint(str_aux2);

            switch( report_format )
                case {'pdf' 'tiff' }

                    if( jj == 1 )
                        export_fig(filename, '-nocrop', '-bookmark', Seconds2HMS(win_length), '-append', ['-' report_format], fig_hdl);
                    else
                        export_fig(filename, '-nocrop', '-bookmark', ['    ' Seconds2HMS(start_time)], '-append', ['-' report_format], fig_hdl);
                    end

                    
                otherwise
                    [aux_path, aux_filename] = fileparts(filename);
                    aux_filename = [aux_path filesep aux_filename];
                    aux_filename = sprintf( '%s_%04d.%s', aux_filename, ii, report_format);
                    export_fig(aux_filename, '-nocrop', ['-' report_format], fig_hdl);
            end

            pb1.checkpoint(str_aux2);

            ii = ii + 1;
            jj = jj + 1;

            pb1.end_loop();
            
        end

        pb.end_loop();

    end

    if( ishandle(fig_hdl) )
        close(fig_hdl)
    end
    
%     clear pb

catch MException
    %% Error handling

    %% No User interface report clearly to log file

    disp_string_framed(2, sprintf('ERROR in rec %s', ECG_w.ECG_header.recname))

    rethrow(MException)

end


function fig_hdl = prepare_fig_hdl( fig_hdl )

    paper_size = [0 0 29.7 21.0 ];

    if( ishandle(fig_hdl) )
        close(fig_hdl)
    end
    fig_hdl = figure('Visible', 'off');
%     fig_hdl = figure();
    set(fig_hdl, 'NumberTitle', 'off');
%     prev_units = get(fig_hdl, 'PaperUnits');
    set(fig_hdl, 'PaperUnits', 'centimeters');
    set(fig_hdl, 'PaperType', 'A4');
    set(fig_hdl, 'PaperOrientation', 'landscape');
    set(fig_hdl, 'PaperPosition', paper_size);
%     disp(get(fig_hdl, 'PaperPosition'));
%     set(fig_hdl, 'PaperUnits', prev_units);
    
%     prev_units = get(fig_hdl, 'Units');
    set(fig_hdl, 'Units', 'centimeters');
    set(fig_hdl, 'Position', paper_size);
%     disp(get(fig_hdl, 'Position'));
%     set(fig_hdl, 'Units', prev_units);
    
    
