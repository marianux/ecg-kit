function [ x_int, t_int, indexes ] = extract_interval( x, t, int_ini, int_end )
%EXTRACT_INTERVAL     Very simple function to extract an interval from a signal
%
% Created by Jesús Lázaro <jlazarop@unizar.es> in 2011
%--------
%   Sintax: [ x_int, t_int, indexes ] = extract_interval( x, t, int_ini, int_end )
%   In:   x = signal
%         t = time vector
%         int_ini = interval begin time (same units as 't')
%         int_end = interval end time (same units as 't')
%
%   Out:  x_int = interval [int_ini, int_end] of 'x'
%         t_int = interval [int_ini, int_end] of 't'
%         indexes = indexes corresponding to returned time interval

    if nargin<4
        error('Not enough input arguments');
    end
    
    indexes = find(t>=int_ini & t <=int_end);
    x_int = x(indexes);
    t_int = t(indexes);
end