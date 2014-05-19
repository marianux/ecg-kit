function UpdStructure = updstruct(OldStructure, NewStructure, SkipEmptyUpdatesOfTheExistingFields)
%
% Structure update
%
%   UpdStructure = updstruct(OldStructure, NewStructure)
%   
% In the result UpdStructure all new filed from NewStructure will be added, 
% coinciding fields will be replaced, fields which do not exist in the NewStruct
% will be preserved

% Copyright: S.Verzakov, s.verzakov@tudelft.nl

if nargin < 3 | isempty(SkipEmptyUpdatesOfTheExistingFields)
  SkipEmptyUpdatesOfTheExistingFields = 0;
end

if nargin < 2 
  NewStructure = [];
end

if ~isempty(NewStructure) & (~isstruct(NewStructure) | length(NewStructure) ~=1)
  error('NewStructure parameter has to be a structure (of length 1) or be empty');
end

if ~isstruct(OldStructure) | length(OldStructure) ~=1
  error('OldStructure parameter has to be a structure (of length 1)');
end

UpdStructure = OldStructure;

if isempty(NewStructure)
  return;
end

fn = fieldnames(NewStructure);
for i=1:length(fn)

  v = getfield(NewStructure,fn{i});

  if ~(SkipEmptyUpdatesOfTheExistingFields & isfield(UpdStructure, fn{i}) & isempty(v))
    if isa(v,'struct') & length(v) == 1
      if isfield(UpdStructure, fn{i}) 
        w = getfield(UpdStructure, fn{i});  
        if isa(w,'struct') & length(w) == 1
          v = updstruct(w,v); 
        end  
      end
    end 
    UpdStructure = setfield(UpdStructure,fn{i},v);
  end
end

return;
