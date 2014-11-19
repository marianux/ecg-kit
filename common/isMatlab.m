function bAux = isMatlab()

bAux = false;

matlab_ver = ver('Matlab');

if( ~isempty(matlab_ver) && strcmpi(matlab_ver.Name, 'MATLAB') )
    bAux = true;
end
