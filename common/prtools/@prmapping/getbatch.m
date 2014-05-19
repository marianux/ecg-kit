%GETBATCH Test on possible execution of mapping in batch mode
%
%    [FLAG,BATCHSIZE,OBJSIZE] = GETBATCH(W)
%
% INPUT
%    W          Mapping
%
% OUTPUT
%    FLAG       Flags batch processing on (TRUE) or off (FALSE)
%    BATCHSIZE  Number of objects per batch
%    OBJSIZE    Number of objects above which batch processing is enabled
%
%DESCRIPTION
%FLAG is TRUE in case batch mode is allowed (default).
%Use SETBATCH to set or disable batch processing of mappings
%
%SEE ALSO
%MAPPINGS, SETBATCH
