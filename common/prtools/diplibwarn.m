%DIPLIBWARN  Warn no DIPimage once

function diplibwarn

persistent DIPWARN

if isempty(DIPWARN)
  DIPWARN = 0;
  disp([newline '--> The DIPIMAGE package is not available. Replacements ' ...
    newline '--> will be used and may behave slightly different.' newline])
end

return