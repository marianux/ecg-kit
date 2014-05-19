% Matlab mex-file
% Function
% 	writeannot
% Purpose
%	Write annotation files for biomedical signals in MIT Format
% Synopsis
%       writedannot(ann_filename,ann_struct)
% Description
% Input ann_struct is a struct with 6 fields: time, anntyp, subtyp, chan, num
% and aux, as the one returned by readannot.  
% Each of this fields is a vector with as many elements as the number of events
% in the anotator.  The struct can be incomplete, but time and anntyp fields are
% strictly necessary.
%
% Copyright (c) Juan Pablo Martinez Cortes. University of Zaragoza
% e-mail: juanpabl@tsc1.cps.unizar.es
