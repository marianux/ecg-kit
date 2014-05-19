%ADDPREPROC Add commands in the PREPROC field of a datafile
%
%   A = ADDPREPROC(A,PREPROC,PARS,OUTSIZE)
%
% INPUT
%   A       - Datafile
%   PREPROC - String with preprocessing command
%   PARS    - Cell array with parameters (default empty)
%   OUTSIZE - Output size of objects
%
% OUTPUT
%   A       - Datafile
%
% DESCRIPTION
% Extends the structure array of preprocessing commands, A.PREPROC( )
% with one element, having two fields:
% A.PREPROC(N).PREPROC = PREPROC
% A.PREPROC(N).PARS = PARS
%
% Preprocessing can only be defined for raw datafiles that do not contain
% datasets. It is an error to define preprocessing for a datafile that
% points to MAT files. In that case datasets are expected and preprocessing
% is superfluous.
%
% The first command in the PREPROC field should always be a file read
% command. The DATAFILE constructor stores by default IMREAD. Use
% SETPREPROC to clear the PREPROC field and replace it by another
% command.
%
% PRTools needs to know the size of the output objects. If not supplied in
% the call, it is found by processing the first object. This may take some
% time.
%
% SEE ALSO 
% DATAFILES, SETPREPROC, ADDPOSTPROC, DATASETS
