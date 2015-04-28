%% (Internal) Write annotation files for biomedical signals in MIT Format. (MEX file)
%
% Description
% Input ann_struct is a struct with 6 fields: time, anntyp, subtyp, chan, num
% and aux, as the one returned by readannot.  
% Each of this fields is a vector with as many elements as the number of events
% in the anotator.  The struct can be incomplete, but time and anntyp fields are
% strictly necessary.
% 
%     writeannot(ann_filename,ann_struct)
% 
% Arguments:
% 
%   + ECG_header: original header struct.
% 
%   + ECG_idx: signal indexes to select.
% 
% Output:
% 
%   + ECG_header: trimed header struct.
% 
% See also read_ecg, writeheader, ECGwrapper
% 
% Author: Juan Pablo Martinez Cortes
% adapted by Mariano Llamedo Soria
% Version: 0.1 beta
% Birthdate: 21/7/2010
% Last update: 20/02/2013
% Copyright 2008-2015
% 

