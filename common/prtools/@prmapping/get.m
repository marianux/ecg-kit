%GET Get mapping parameter fields
%
%   [VALUE1,VALUE2,...] = GET(W,FIELD1,FIELD2,...)
%
% INPUT
%   W       Mapping
%   FIELDx  Field names (strings)
% 
% OUTPUT
%   VALUEx  Value of field x
%
% DESCRIPTION
% Get the values of the fields in W specified by FIELD1 etc. Note that mapping
% fields may also be extracted directly by, e.g. W.DATA.
%
% List of field names (see also MAPPING):
%
%   MAPPING_FILE  name of the routine used for learning or executing the mapping
%   MAPPING_TYPE  string defining the type of mapping:
%                 'untrained', 'trained', "combiner' or 'fixed'.
%   DATA          data, structure or cell array needed for defining the mapping.
%   LABELS        array with labels to be used as feature labels for the dataset
%                 created by executing the mapping.
%   SIZE_IN       input dimensionality or size vector describing its shape.
%   SIZE_OUT      output dimensionality or size vector describing its shape.
%   SCALE         output multiplication factor or vector.
%   OUT_CONV      0,1,2,3 for defining the desired output conversion:
%                 0 - no (default), 1 - SIGM, 2 NORMM or 3 - SIGM and NORMM.
%   COST          classification cost matrix
%   NAME          string with mapping name.
%   USER          user definable variable.
%   VERSION       version field
%
% EXAMPLES
% [DATA,LABELS] = GET(W,'data','labels')
%
% SEE ALSO
% DATASETS, MAPPINGS, MAPPING
