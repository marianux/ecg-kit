function str_aux = disp_option_enumeration(srtMessage, cOptionEnum)
% This function format a string for fprintf like functions in order to
% present a message srtMessage followed by an enumeration of string
% elements included in the cell str cOptionEnum

str_aux = [ repmat(' + ', length(cOptionEnum), 1) char(cOptionEnum) repmat('\n', length(cOptionEnum), 1 ) ];
str_aux = [ srtMessage '\n' rowvec(str_aux') ];
