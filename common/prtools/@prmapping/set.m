%SET Set mapping parameters
%
%	    W = SET(W,NAME1,VALUE1,NAME1,VALUE1,...)
%
% Sets the fields given by their names (as strings) of the mapping W to the
% supplied values.E.G.: W = SET(W,'data',{DATA1,DATA2},'labels',LABELS).
% This is not different from using the field specific routines
% (e.g. SETDATA(W,{DATA1,DATA2})
%
% List of field names, see also MAPPING:
%
% MAPPING_FILE  name of the routine used for learning or executing the mapping
% MAPPING_TYPE  string defining the type of mapping:
%               'untrained', 'trained', "combiner' or 'fixed'.
% DATA          Data, structure or cell array needed for defining the mapping.
% LABELS        Array with labels to be used as feature labels for the dataset
%               created by executing the mapping.
% SIZE_IN       Input dimensionality or size vector describing its shape.
% SIZE_OUT      Output dimensionality or size vector describing its shape.
% SIZE          Not a field, but sets both, SIZE_IN and SIZE_OUT from a vector.
% SCALE         Output multiplication factor or vector.
% OUT_CONV      0,1,2,3 for defining the desired output conversion:
%               0 - no (default), 1 - SIGM, 2 NORMM or 3 - SIGM and NORMM.
% COST          classification cost matrix.
% NAME          String with mapping name.
% USER          User definable variable.
%
% See also DATASETS, MAPPINGS, MAPPING
