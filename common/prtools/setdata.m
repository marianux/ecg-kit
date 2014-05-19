%SETDATA Dummy routine for setdata for doubles
%
%   OUT = SETDATA(A,B,PARS)
%
%DESCRIPTION
%This routine catches calls to SETDATA in case A is not a dataset but a
%double. In that case just B will be returned instead of a dataset with all
%annotations as defined by A but with data B.

function b = setdata(a,b,varargin)

return