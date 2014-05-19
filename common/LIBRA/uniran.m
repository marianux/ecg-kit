function [random,seed]=uniran(seed)

%UNIRAN is the uniform random generator used in mcdcov.m and ltsregres.m
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html

seed=floor(seed*5761)+999;
quot=floor(seed/65536);
seed=floor(seed)-floor(quot*65536);
random=seed/65536.D0;
