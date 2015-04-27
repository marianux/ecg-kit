%% (Internal) Display a message framed to a string
%   
%     This function format a string for fprintf like functions in order to
%     present a message string2display framed with # characters
% 
%   str_out = disp_string_framed(style, string2display )
% 
% Arguments:
% 
%      + style: if numeric, could be the stream identifier for fprintf
%      function. If char, can be the "style" argument for cprintf function.
% 
%      + string2display: The string to frame. 
% 
% Output:
% 
%      + str_out: resulting string
% 
% Example:
% 
% 
% See also disp_option_enumeration, disp_string_title
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
function str_out = disp_string_framed(style, string2display )

    str_out = [];
    
    if( nargin == 1 || isempty(string2display) )
        string2display = style;
        style = 'text';
        fid = 1;
    end

    if( isempty(style) )
        style = 'text';
        fid = 1;
    elseif( isnumeric(style) )
        fid = style;
        style = 'text';
    else
        if( ischar(style) )
            fid = 1;
        else
            style = 'text';
            fid = 1;
        end
    end

    str_aux = sprintf('# %s #', string2display);
    str_aux1 = repmat('#', 1, length(str_aux));
    if( fid == 0)
        string2display = strrep(string2display, '\', '\\');
        str_out = sprintf('\n%s\n# %s #\n%s\n\n', str_aux1, string2display, str_aux1);
    elseif( fid == 1)
        fprintf(fid, '\n%s\n# ', str_aux1);
        string2display = strrep(string2display, '\', '\\');
        cprintf(style, [string2display ' ']);
        fprintf(fid, '#\n%s\n\n', str_aux1);
    else
        fprintf(fid, '\n%s\n# %s #\n%s\n\n', str_aux1, string2display, str_aux1);
    end

end
