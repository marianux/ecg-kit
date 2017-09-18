%% (Internal) function reads the header of signal files
%
%     heasig = readheader(hea_full_filename)
% 
% Arguments:
% 
%	Input parameters: character string with hea_full_filename of header file
% 
% 
% Output:
% 
%	Output parameter: struct heasig with header information
% 
% 
% Example:
% 
% 
% See also read_ECG
% 
% Authors:
% Salvador Olmos
% Last update Juan Bolea 07/10/2008
% adapted to the kit by Mariano Llamedo Soria (llamedom at frba.utn.edu.ar)
% Version: 0.1 beta
% Birthdate: 29/05/2014
% Last update: 29/05/2014
% Copyright 2008-2015
%
function heasig = readheader(hea_full_filename)


heasig = [];

try 
    fid = fopen(hea_full_filename,'rt');
    if (fid <= 0) % Rute 11/08/2010
        fid = fopen([hea_full_filename '.hea'],'rt');%Rute 11/08/2010
    end
    if (fid <= 0)
        %     heasig = [];
        return;
        %     disp(['error in opening file ' hea_filename]);
    end
    [hea_path hea_filename hea_ext] = fileparts(hea_full_filename);

    files = dir([hea_path filesep hea_filename '.*'] );
    NotHea_idx = find(~strcmpi( [hea_filename hea_ext], {files(:).name}));
    [~, max_idx] = max(cell2mat({files(NotHea_idx).bytes}));

    if( isempty(max_idx) )
        warning('readheader:NoECGfound', ['Could not found any ECG signal in the same path of ' hea_full_filename ])
        ecg_filename = 'ECG_file_not_found';
    else
        ecg_filename = files(NotHea_idx(max_idx)).name;
    end

    % Used symbols to distinguish some fields of the header file
    pp = ' /+:()x';


    % First line reading
    s = fgetl(fid);
    % Remove blank or commented lines
    while s(1) == '#'
        s = fgetl(fid);
    end

    %% Record line

    % Example:  hea_filename'/'number_of_segments[opt] number_of_signals
    % sampling_frequency[opt]'/'counter_frequency[opt]'('base_counter_value')'[opt]
    % number_of_samples_per_signal[opt] base_time[opt] base_date[opt]

    %%
    % record hea_filename
    [heasig.recname,s] = strtok(s,pp);
    % multi-segment record
    if strcmp(s(1),'/')
        [s1 s] = strtok(s,pp);
        heasig.nsegm = str2double(s1);
    end
    % number of signals
    [s1 s] = strtok(s,pp);
    heasig.nsig = str2double(s1);

    if(heasig.nsig < 1 || heasig.nsig > 20 )

        warning('readheader:AtypicalNsig', ['Atypical number of signals ' num2str(heasig.nsig) ])

        if( heasig.nsig == 0 )
            % NSRDB have nsig = 0
            heasig.nsig = 1;
        end
    end

    % sampling frequency (in samples per second per signal) [optional]
    if ~isempty(s)
        [s1 s] = strtok(s,pp);
        heasig.freq = str2double(s1);
        if heasig.freq > 0 && heasig.freq < 1
            heasig.freq = 1/heasig.freq;
        end
    end

    if isfield(heasig,'freq')  % if exists sampling frequency field
        if ~isempty(deblank(s))
            if strcmp(s(1),'/')
                [s1 s] = strtok(s,pp);  % Counter frequency (in ticks per second) [optional]
                heasig.cntfreq = str2double(s1);
                if strcmp(s(1),'(')   % (453) 765 834
                    [s1 s] = strtok(s,'('); % 453) 765 834
                    [s1 s] = strtok(s1,')'); % s1 = 453 s = ) 765 834
                    heasig.basecnt = str2double(s1);  % Base counter value [optional]
                    [s1 s] = strtok(s,pp); % remove parentheses
                    heasig.nsamp = str2double(s1);
                else
                    [s1 s] = strtok(s,pp);  %20ABR2010 Juan
                    heasig.nsamp = str2double(s1); % number of samples per signal [optional]
                end
            else
                [s1 s] = strtok(s,pp);
                heasig.nsamp = str2double(s1); % number of samples per signal [optional]

                if( heasig.nsamp == 0 )
                    % NSRDB have nsamp = 0
                    warning('readheader:AtypicalNsig', ['Atypical number of samples ' num2str(heasig.nsamp) ])
                end

            end
        end
        if isfield(heasig,'nsamp')  % if exists number of samples per signal
            if isempty(deblank(s)) || isempty(strfind(s,':'))
                heasig.btime = '00:00:00';
                heasig.bdate = '01/01/2000';
            else            
                [s1 s] = strtok(s,':');
                hour = strtrim(s1);
                [s1 s] = strtok(s,':');
                min = strtrim(s1);
                [s1 s] = strtok(s,pp);
                sec = strtrim(s1);
                heasig.btime = [hour ':' min ':' sec];
                if isempty(deblank(s)) || isempty(strfind(s,'/'))
                    heasig.bdate = '01/01/2000';
                else
                    [s1 s] = strtok(s,'/');
                    day = strtrim(s1);
                    [s1 s] = strtok(s,'/');
                    month = strtrim(s1);
                    s1 = strtok(s,pp);
                    year = strtrim(s1);
    %                 yy = datevec(['01/01/' year]);
                    heasig.bdate = [day '/' month '/' year];
                end
            end
        end
    end

    if isfield(heasig,'nsegm')
        %% Segment Specification Lines    
        %Example:  record_name number_of_samples_per_signal
        j = 1;
    %     addt = 0;
        for i = 1:heasig.nsegm
            s = fgetl(fid);
            [segname{i},s] = strtok(s,pp);
            if i > 1
                heasig.sumsampsegm(i) = str2double(s)+heasig.sumsampsegm(i-1);
            else
                heasig.sumsampsegm(i) = str2double(s);
            end
            heasig.sampsegm(i) = str2double(s);
            if ~strcmp(segname{i}(1),'~')
                hsig{i} = readheader([hea_path filesep segname{i}]);
                for k = 1:hsig{i}.nsig
                    desc{j} = strtrim(hsig{i}.desc(k,:));
                    j = j+1;
                end
            else
                hsig{i}.null = 1;
                hsig{i}.nsamp = heasig.sampsegm(i);
            end
        end    
        heasig.fname = char(segname); 
        heasig.desc = char(unique(desc));
        heasig.hd = hsig;
    %     clear heasig;
    %     heasig = hsig;        
    else

        %% Signal Specification Lines

        %Example:  file_hea_filename format'x'samples_per_frame[opt]':'skew[opt]'+'byte_offset[opt]
        % ADC_gain[opt]'('baseline_ADC_units')'[opt]'/'units[opt]
        % ADC_resolution[opt] ADC_zero[opt] initial_value[opt] checksum[opt]
        % block_size[opt] description[opt]

        %%
        heasig.spf = ones(1,heasig.nsig);
        heasig.baseline = zeros(1,heasig.nsig);
        heasig.units(1:heasig.nsig) = {''};
        heasig.fname = cell(heasig.nsig,1);
        heasig.group = nan(heasig.nsig,1);
        heasig.fmt = nan(heasig.nsig,1);
        heasig.offset = nan(heasig.nsig,1);
        heasig.gain = nan(heasig.nsig,1);
        heasig.adcres = nan(heasig.nsig,1);
        heasig.adczero = nan(heasig.nsig,1);
        heasig.initval = nan(heasig.nsig,1);
        heasig.baseline = nan(heasig.nsig,1);
        heasig.adczero = nan(heasig.nsig,1);
        heasig.cksum = nan(heasig.nsig,1);
        heasig.bsize = nan(heasig.nsig,1);
        heasig.spf = nan(heasig.nsig,1);
        heasig.desc = num2str((1:heasig.nsig)');

        % Reading nsig lines, corresponding one for every lead
        for i = 1:heasig.nsig
            if( feof(fid) )
                break
            else
                s = fgetl(fid);
            end

            % Remove blank or commented lines
            while s(1)=='#'
                if( feof(fid) )
                    break
                else
                    s = fgetl(fid);
                end
            end

            if( s(1)=='#' && feof(fid) ); break; end;

            % file hea_filename
            [heasig.fname{i},s] = strtok(s);  %% Modificado 29/04/2008  old "[heasig.fname(i,:),s]=strtok(s,pp);"
            [s1,s] = strtok(s,pp);

            %     % group
            if i == 1
                heasig.group(i) = 0;
            else
                if strcmp(heasig.fname(i,:),heasig.fname(i-1,:))
                    heasig.group(i)=0;
                else
                    heasig.group(i) = heasig.group(i-1) + 1;
                end
            end
            % format
            [s3,s4]=strtok(s);
            a=[strfind(s3,'x') strfind(s3,':') strfind(s3,'+')];
            if isempty(a)
                heasig.fmt(i)=str2double(deblank(s1));
            else
                [s2,s]=strtok(s);
                a=[a length(s2)+1];
                for k=1:length(a)-1
                    switch (s2(a(k)))
                        case 'x' % samples per frame [optional]
                            heasig.fmt(i)=str2double(s1);
                            heasig.spf(i)=str2double(s2(a(k)+1:a(k+1)-1));
                        case ':' % skew [optional]
                            heasig.fmt(i)=str2double(s1);
                            heasig.skew(i)=str2double(s2(a(k)+1:a(k+1)-1));
                        case '+' % byte offset [optional]
                            heasig.fmt(i)=str2double(s1);
                            heasig.offset(i)=str2double(s2(a(k)+1:a(k+1)-1));
                    end
                end
            end
            [s1,s]=strtok(s); %% Modificado 29/04/2008  old "[s1,s]=strtok(s,pp);"
            if ~isempty(deblank(s))
                a=[strfind(s1,'(') strfind(s1,'/')];

                %%  Modificado 29/04/2008
                if isempty(a) % ADC gain (ADC units per physical unit) [optional] -> units assumed mV
                    heasig.gain(i) = str2double(s1);
                    heasig.units{i} = 'mV';
                else
                    [s2,s1] = strtok(s1,pp);
                    heasig.gain(i) = str2double(s2);
                    if ~isempty(s1)
                        a = strfind(s1,'/');
                        if isempty(a)  % baseline (ADC units) [optional]
                            heasig.baseline(i) = str2double(s1(2:end-1));
                            %                 heasig.units(i,:) = 'mV';
                        else
                            if strcmp(s1(1),'(')
                                heasig.baseline(i) = str2double(s1(2:a-2));
                                f = s1(a+1:end);
                                heasig.units{i} = f; % units[optional]
                            else
                                heasig.baseline(i) = 0;
                                f = s1(2:end);
                                heasig.units{i} = f;
                            end
                        end
                    end
                end
            end

            if isnan(heasig.gain(i))
                heasig.gain(i) = 200;
            else
                if heasig.gain(i) == 0
                    heasig.gain(i) = 200;
                end
                if ~isempty(deblank(s))
                    s_aux = s;
                    [s1,s] = strtok(s,pp);
                    heasig.adcres(i) = str2double(s1); % ADC resolution (bits) [optional]
                    if isnan(heasig.adcres(i))
                        % failed to parse this field
                        s = s_aux ;
                        heasig.adcres(i) = 12;
                    elseif (heasig.adcres(i) == 0)
                        heasig.adcres(i) = 12;
                    end

                end
                if ~isnan(heasig.adcres(i))
                    if ~isempty(deblank(s))
                        s_aux = s;
                        [s1,s] = strtok(s,pp);
                        heasig.adczero(i) = str2double(s1); % initial value [optional]
                        if isnan(heasig.adczero(i))
                            % failed to parse this field
                            s = s_aux ;
                        end
                    end
                end

                if ~isnan(heasig.adczero(i))
                    if ~isempty(deblank(s))
                        s_aux = s;
                        [s1,s] = strtok(s,pp);
                        heasig.initval(i) = str2double(s1); % initial value [optional]
                        if isnan(heasig.initval(i))
                            % failed to parse this field
                            s = s_aux ;
                        end
                    end
                end
                if ~isnan(heasig.initval(i))
                    if ~isempty(deblank(s))
                        s_aux = s;
                        [s1,s] = strtok(s,pp);
                        heasig.cksum(i) = str2double(s1); % initial value [optional]
                        if isnan(heasig.cksum(i))
                            % failed to parse this field
                            s = s_aux ;
                        end

                    end
                end
                if ~isnan(heasig.cksum(i))
                    if ~isempty(deblank(s))
                        s_aux = s;
                        [s1,s] = strtok(s,pp);
                        heasig.bsize(i) = str2double(s1); % initial value [optional]
                        if isnan(heasig.bsize(i))
                            % failed to parse this field
                            s = s_aux ;
                        end
                    end
                end
                if  ~isnan(heasig.bsize(i))
                    if ~isempty(strtrim(s))
                        heasig.desc(i,1:length(s))=s; % description [optional]
                    end
                end
            end

        end
        if isfield(heasig,'units')
            heasig.units = char(heasig.units);
        end

        heasig.fname = char(heasig.fname);

    end
    fclose(fid);

catch ME
    fclose(fid);
end
