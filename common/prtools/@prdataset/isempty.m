%ISEMPTY Dataset overload
%
%	I = ISEMPTY(A,FIELD)
%
% INPUT
%  A     Dataset
%  FIELD Dataset field, default 'data'
%
% OUTPUT
%  I     Flag, 1 if field is empty, 0 otherwise. 
%
% DESCRIPTION
% Dataset overload for ISEMPTY. This is particulary useful for
% ISEMPTY(A) to test on an empty dataset, and
% ISEMPTY(A,'prior') to test on an undefined PRIOR field
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% PRDATASET, SETPRIOR, GETPRIOR
