%% (Internal) Create a string with the invocation of a function
%   
% Example
% 
%   strAux = GetFunctionInvocation( funcName, args)
% 
% Arguments:
% 
%      + funcName: Confusion matrix of n_classes x n_classes
% 
%      + args: a cell with the arguments used 
% 
% Output:
% 
%     + strAux the string with the function invocation.
% 
% Example:
% 
%       strAux = GetFunctionInvocation(mfilename, varargin);
% 
% See also a2hbc_main
% 
% Author: Mariano Llamedo Soria llamedom@electron.frba.utn.edu.ar
% Version: 0.1 beta
% Last update: 14/5/2014
% Birthdate  : 21/4/2015
% Copyright 2008-2015

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
