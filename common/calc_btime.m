function btime = calc_btime( start_sample, sampling_rate )

hours = floor((start_sample-1)/sampling_rate/60/60);
mins = floor((start_sample-1)/sampling_rate/60 - hours *60);
secs = floor((start_sample-1)/sampling_rate - mins *60 - hours * 60 * 60);
milli = round(((start_sample-1)/sampling_rate - mins *60 - hours * 60 * 60 - secs) * 1000);
btime = sprintf('%0d:%0d:%0d:%03d',hours, mins, secs, milli);

