function [vec,I]=greatsort(x);

%GREATSORT sorts the vector x in descending order.
%
% Required input arguments: 
%     x : vector to be sorted
%
% Output arguments:
%   vec : sorted vector
%    I  : index => x(I)=vec.
%
% I/O: [vec,I]=greatsort(x);
% 
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html

vec=-x;
[svec,I]=sort(vec);
vec=-svec;

