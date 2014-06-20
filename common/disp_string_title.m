function disp_string_title(stream_id, string2display )

    str_aux = sprintf(' %s ', string2display);
    str_aux1 = repmat('-', 1, length(str_aux));
    fprintf(stream_id, '\n%s\n%s\n\n', str_aux, str_aux1);

end