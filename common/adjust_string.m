function adjusted_str = adjust_string(str2trimm, target_width, where2trimm )

    if( nargin < 3 || isempty(where2trimm) )
       where2trimm = 'center';
    end

    if( nargin < 2 || isempty(target_width) || target_width < 7 )
       target_width = 20;
    end

    missing_chars_symbol = ' ... ';
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
        adjusted_str = [str2trimm repmat(' ', 1, target_width - lstr2trimm)];
    end
