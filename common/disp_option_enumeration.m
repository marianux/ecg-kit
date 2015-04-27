%% (Internal) Display an enumeration of options to a string
%   
%     This function format a string for fprintf like functions in order to
%     present a message srtMessage followed by an enumeration of string
%     elements included in the cell str cOptionEnum
% 
%   str_aux = disp_option_enumeration(srtMessage, cOptionEnum)
% 
% Arguments:
% 
%      + srtMessage: A title message previous to the enumeration
% 
%      + cOptionEnum: cellstring with the items to enumerate 
% 
% Output:
% 
%      + str_aux: resulting string
% 
% Example:
% 
% 
% See also disp_string_framed, disp_string_title
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
function str_aux = disp_option_enumeration(srtMessage, cOptionEnum)

%  Escape problematic strings, such as file paths 
cOptionEnum = cellfun(@(a)(strrep(a, '\', '\\')), cOptionEnum, 'UniformOutput', false );

str_aux = [ repmat(' + ', length(cOptionEnum), 1) char(cOptionEnum) repmat('\n', length(cOptionEnum), 1 ) ];
str_aux = [ srtMessage '\n' rowvec(str_aux') ];
