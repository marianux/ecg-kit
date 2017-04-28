%% (Internal) Works with strings to center, trim and justify to a certain string width
%   
%   adjusted_str = adjust_string(str2trimm, target_width, where2trimm )
% 
% Arguments:
% 
%      + str2trimm: the string
% 
%      + target_width: size of the target string.
% 
%      + where2trimm: "left" "right" "center"
% 
% Output:
% 
%      + adjusted_str: resulting string.
% 
% Example:
% 
% See also addpath
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015
% 
function adjusted_str = adjust_string(str2trimm, target_width, where2trimm, missing_chars_symbol )

    if( nargin < 4 || isempty(missing_chars_symbol) )
        missing_chars_symbol = ' ... ';
    end

    if( nargin < 3 || isempty(where2trimm) )
       where2trimm = 'center';
    end

%     if( nargin < 2 || isempty(target_width) || target_width < 7 )
    if( nargin < 2 || isempty(target_width) )
       target_width = 20;
    end

    lmissing_chars_symbol =length(missing_chars_symbol);
    
    lstr2trimm = length(str2trimm);

    if( lstr2trimm > target_width )
        
        sample_length = target_width - lmissing_chars_symbol;
        
        switch(lower(where2trimm))

            case 'left'
                adjusted_str = [ str2trimm(end-sample_length+1:end) missing_chars_symbol] ;

            case 'center'
                left_length = sample_length - round(sample_length/2);
                right_length = sample_length - left_length;
                adjusted_str = [str2trimm(1:left_length) missing_chars_symbol str2trimm(end-right_length+1:end) ];

            case 'right'
                adjusted_str = [missing_chars_symbol str2trimm(1:sample_length)];

            otherwise
                error('Dont understand "casewhere2trimm". "left" "right" "center" ')
        end
    else
        
        switch(lower(where2trimm))
            case 'pad'
                adjusted_str = [str2trimm repmat(' ', 1, target_width - lstr2trimm)];
                
            case 'center'
                aux_val = floor((target_width - lstr2trimm)/2);
                adjusted_str = [repmat(' ', 1, aux_val ) str2trimm repmat(' ', 1, target_width - lstr2trimm - aux_val )];
                
            case 'right'
                adjusted_str = [repmat(' ', 1, target_width - lstr2trimm) str2trimm ];
                
            case 'left'
                adjusted_str = [str2trimm repmat(' ', 1, target_width - lstr2trimm) ];
                
            otherwise
                adjusted_str = str2trimm;
        end
        
    end
