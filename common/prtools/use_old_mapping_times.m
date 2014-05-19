%USE_OLD_MAPPING_TIMES
%
%    USE_OLD_MAPPING_TIMES(FLAG)
%
% INPUT
%   FLAG  TRUE or FALSE
%
% DESCRIPTION
% For FLAG = TRUE (1) the old version of the multiplication of mappings
% (e.g. W = 100*FISHERC) will be used. This may be reset by a call with
% FLAG = FALSE

function use_old_mapping_times(flag)

global OLD_MAPPING_TIMES;
if flag
  OLD_MAPPING_TIMES = true;
else
  OLD_MAPPING_TIMES = false;
end