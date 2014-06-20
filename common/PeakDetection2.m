function peakloc = PeakDetection2(data,fs,varargin)
%
% peaks = PeakDetection2(x,fs,wlen,fp1,fp2,th,flag),
% R-peak detector based on Pan-Tompkins method.
%
% inputs:
% x: vector of input data
% fs: sampling rate
% wlen: moving average window length (default = 150ms)
% fp1: lower cut-off frequency (default = 10Hz)
% fp2: upper cut-off frequency (default = 33.3Hz)
% th: detection threshold (default = 0.2)
% flag: search for positive (flag=1) or negative (flag=0) peaks. By default
% the maximum absolute value of the signal, determines the peak sign.
%
% output:
% peaks: vector of R-peak impulse train
%
%
% Open Source ECG Toolbox, version 2.0, March 2008
% Released under the GNU General Public License
% Copyright (C) 2008  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- GIPSA-Lab, INPG, Grenoble, France
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

if(nargin>2  && ~isempty(varargin{1})),
    winlen = varargin{1};
else
    winlen = .150; % 150ms
end

if(nargin>3 && ~isempty(varargin{2})),
    fp1 = varargin{2};
else
    fp1 = 10;
end

if(nargin>4  && ~isempty(varargin{3})),
    fp2 = varargin{3};
else
    fp2 = 33.3;
end

if(nargin>5  && ~isempty(varargin{4})),
    thr = varargin{4};
else
    thr = 0.2;
end

if(nargin>6  && ~isempty(varargin{5})),
    flag = varargin{5};
else
    flag = abs(max(data))>abs(min(data));
end


N = length(data);
data = data(:);

L1 = round(fs/fp2);    % first zero of the LP filter is placed at f = 33.3Hz;
L2 = round(fs/fp1);    % first zero of the HP filter is placed at f = 3Hz;

fprintf(1, '.');

% x0 = data - Median(data',N,round(fs*winlen/3));

q = ceil(fs/150);
N_padded = (ceil(N / q) * q);
data = [data; zeros(N_padded-N,1)];
data = data - resample(MedianFilt(resample(data,1,q), round(fs*winlen/3/q)), q, 1);
data = data(1:N);
% x0 = data;

fprintf(1, '.');

% LP filter
x = filter([1 zeros(1,L1-1) -1],[L1 -L1],data);
x = filter([1 zeros(1,L1-1) -1],[L1 -L1],x);
x = [x(L1:end);zeros(L1-1,1) + x(end)]; % lag compensation

fprintf(1, '.');

% HP filter
y = filter([L2-1 -L2 zeros(1,L2-2) 1],[L2 -L2],x);

% differentiation
z = diff([y(1) ; y]);

% squaring
w = z.^2;

fprintf(1, '.');

% moving average
L3 = round(fs*winlen);
v = filter([1 zeros(1,L3-1) -1],[L3 -L3],w);
v = [v(round(L3/2):end);zeros(round(L3/2)-1,1) + v(end)]; % lag compensation

vmax = max(v);
p = v > (thr*vmax);

fprintf(1, '.');

% edge detection
rising  = find(diff([0 ; p])==1);      % rising edges
falling = find(diff([p ; 0])==-1);     % falling edges

if( length(rising) == length(falling)-1 ),
    rising = [1 ; rising];
elseif( length(rising) == length(falling)+1 ),
    falling = [falling ; N];
end

fprintf(1, '.');

peakloc = zeros(length(rising),1);
width = zeros(length(rising),1);

Nrising = max(1, length(rising)/20);
j = 0;

if(flag)
    for i=1:length(rising)
        
        if( j > Nrising)
            fprintf(1, '.');
            j = 0;
        else
            j = j + 1;
        end
        
        [val mx] = max( data(rising(i):falling(i)) );
        peakloc(i) = mx - 1 + rising(i);
        width(i) = falling(i) - rising(i);
    end
else
    for i=1:length(rising)
        if( j > Nrising)
            fprintf(1, '.');
            j = 0;
        else
            j = j + 1;
        end
        
        [val mn] = min( data(rising(i):falling(i)) );
        peakloc(i) = mn - 1 + rising(i);
        width(i) = falling(i) - rising(i);
    end
end
fprintf(1, '\n');

% peaks = zeros(1,N);
% peaks(peakloc) = 1;
