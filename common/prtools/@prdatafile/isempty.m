%ISEMPTY Datafile overload
%
%	I = ISEMPTY(A,FIELD)
%
% INPUT
%  A     Datafile
%  FIELD Datafile field
%
% OUTPUT
%  I     Flag, 1 if field is empty, 0 otherwise. 
%
% DESCRIPTION
% Dataset overload for ISEMPTY. This is particulary useful for
% ISEMPTY(A) to test on an empty datafile, and
% ISEMPTY(A,'prior') to test on an undefined PRIOR field.
%
% A datafile is empty if no files are found in the directory DIR
% after a PRDATAFILE(DIR) definition.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATAFILE
