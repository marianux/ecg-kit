function strAux = GetFunctionInvocation( funcName, args)

strAux = sprintf( [ 'Invocation: ' funcName '(']);
for ii = 1:length(args)
    if( ischar(args{ii}) )
        strAux = [strAux sprintf('''%s''', args{ii})];
    elseif(isempty(args{ii}))
        strAux = [strAux '[]'];
    elseif(isnumeric(args{ii}))
        if( max(size(args{ii})) < 10 )
            strAux = [strAux sprintf( num2str(args{ii}))];
        else
            strAux = [strAux '#BIG_NUMERIC_DATA#'];
        end
    elseif(islogical(args{ii}))
        if( args{ii} )
            strAux = [strAux sprintf('true')];
        else
            strAux = [strAux sprintf('false')];
        end
    end
    strAux = [strAux sprintf(', ')];
end
strAux = [strAux(1:end-2) ')' ];
