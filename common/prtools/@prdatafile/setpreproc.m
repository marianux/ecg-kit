%SETPREPROC Set the PREPROC field of a datafile
%
%   A = SETPREPROC(A,PREPROC,PARS)
%   A = SETPREPROC(A,PSTRUCT)
%   A = SETPREPROC(A)
%
% INPUT
%   A       - Datafile
%   PREPROC - String with preprocessing command
%   PARS    - Cell array with parameters (default empty)
%   PSTRUCT - Structure array with a set of preprocessing commands
%
% OUTPUT
%   A       - Datafile
%
% DESCRIPTION
% Resets the structure array of preprocessing commands in A.PREPROC.
% A.PREPROC(N).PREPROC = PREPROC
% A.PREPROC(N).PARS = PARS
%
% A call without PREPROC and PARS clears A.PREPROC.
%
% Preprocessing can only be defined for raw datafiles that do not contain
% datasets. It is an error to define preprocessing for a datafile that
% points to MAT files. In that case datasets are expected and preprocessing
% is superfluous. All preprocessing commands are executed just before a 
% DATAFILE is converted into a PRDATASET.
%
% The first command in the PREPROC field should always be a file read
% command. The DATAFILE constructor stores by default IMREAD. It is removed
% by SETPREPROC. Be sure to start a new series of preprocessing commands by
% a command to read files. The first parameter of this commands should be
% the filename.
%
% Additional preprocessing commands may be stored by ADDPREPROCC.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>) 
% DATAFILES, ADDPREPROC, DATASETS
